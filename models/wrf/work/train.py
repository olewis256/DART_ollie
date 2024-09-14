import sys
import os
import torch
import random
import numpy as np
import argparse
from torch.utils.data import DataLoader, Dataset
from torch.utils.data.distributed import DistributedSampler
import glob
from torch import nn
import torch.distributed as dist
import iris
import logging
#from networks.afnonet import AFNONet
#from utils.YParams import YParams
from YParams import YParams
import os
# from utils.darcy_loss import LpLoss
from darcy_loss import LpLoss
import time
from torch.optim.lr_scheduler import CosineAnnealingLR
from torch.nn.parallel import DistributedDataParallel
import matplotlib.pyplot as plt
import torch.distributed.elastic.multiprocessing.errors as errors

#import ai_models_fourcastnetv2.fourcastnetv2 as nvs
from modulus.sfnonet import SphericalFourierNeuralOperatorNet as SFNO

def count_parameters():
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def save_checkpoint(model):
 
    torch.save({'model_state': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict()}, 'sfno_ckpt_amp_3.tar')

def restore_checkpoint(checkpoint_path, model):

    checkpoint = torch.load(checkpoint_path)

    model.load_state_dict(checkpoint['model_state'])


    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])


class GetDataset(Dataset):
    def __init__(self, location, RR):
        
        self.location = location
        self.dt = 4
        self.n_history = 1
        self.RR = RR
        self._get_files_stats()
        self.in_chan = 13
        self.out_chan = 12

        self.normalise = True #by default turn on normalization if not specified in config

    def _get_files_stats(self):
        self.files_paths = glob.glob(self.location + "/*.pc*_*") #collect .pc files within omega directory
        self.files_paths.sort()
        self.n_years = len(self.files_paths)
        print("files: ", len(self.files_paths))
        
        self.n_samples_per_cube = 160
        self.img_shape_x = 144
        self.img_shape_y =90
 
        self.n_samples_total = self.n_samples_per_cube * self.n_years

        self.files = [None for _ in range(self.n_years)]

        self.rr = np.full((1,self.img_shape_y, self.img_shape_x), self.RR)
        
    def _open_file(self, cube_idx):
          
        
        data = iris.load(self.files_paths[cube_idx])
        self.files[ cube_idx] = self.files_paths[cube_idx]
        
        regridder = iris.analysis.Linear().regridder(data[0], data[2]) # map (91, 144) cube to (90, 144) to fit with patch size
        
        data[0] = regridder(data[0])
        data[1] = regridder(data[1])
        data[5] = regridder(data[5])
        data[6] = regridder(data[6])
        
#         data[0] = data[0].interpolate(self.points, iris.analysis.Linear())
        
        self.cubes = []
        self.surf = []
        for i in [0,5,6]:
            self.cubes.append(data[i].data[:,4:8,:,:])
#         for i in [2,3]:
#             self.surf.append(data[i].data)  # normalise pressure
        
#         self.cubes[1] /= 1e3
#         self.surf[0] /= 1e2
        
        self.cubes = np.stack(self.cubes)
#         self.surf = np.stack(self.surf)
        
        self.cubes = np.transpose(self.cubes, (1, 0, 2, 3, 4))
        self.cubes = np.reshape(self.cubes, (self.n_samples_per_cube, self.out_chan, self.img_shape_y, self.img_shape_x))
#         self.surf = np.reshape(self.surf, (self.n_samples_per_cube, 2, self.img_shape_y, self.img_shape_x))

#         self.cubes = np.concatenate((self.cubes, self.surf), axis=1)
        
    
        self.output = np.reshape(self.cubes, (self.n_samples_per_cube,self.out_chan, self.img_shape_y, self.img_shape_x)) #quick test, only use 950 hPa temperature

        
    def __len__(self):
        return self.n_samples_total

    def __getitem__(self, global_idx):
 
        cube_idx = int(global_idx/self.n_samples_per_cube) #which 40 day sample file we are in 
        local_idx = int(global_idx%self.n_samples_per_cube) #which 6 hour time we are in within the 40 day sample file
        
        if self.files[cube_idx] is None:
            
            self._open_file(cube_idx)

        y_roll = np.random.randint(0, self.img_shape_x)
      
        sample = self.output[local_idx, :, :, :]
        sample = np.concatenate((sample, self.rr), axis=0)

        if params.two_step_training:
            if local_idx >= self.n_samples_per_cube - 2*self.dt:
                #set local_idx to last possible sample in a cube that allows taking two steps forward
                local_idx = self.n_samples_per_cube - 3*self.dt
                target = self.output[local_idx+self.dt:local_idx+self.dt+2, :, :, :]
                sample = np.reshape(sample, (self.in_chan, self.img_shape_y, self.img_shape_x))
            
            target = self.output[local_idx+self.dt:local_idx+self.dt+2, :, :, :]
            target = np.reshape(target, (self.out_chan, self.img_shape_y, self.img_shape_x))
    
        
        else:
            if local_idx >= self.n_samples_per_cube - self.dt:
            #set local_idx to last possible sample in a cube that allows taking two steps forward
                local_idx = self.n_samples_per_cube - 2*self.dt

            target = self.output[local_idx+self.dt, :, :, :]
    
            target = np.reshape(target, (self.out_chan, self.img_shape_y, self.img_shape_x))
            
        sample = np.roll(sample, y_roll, axis = -1)
        target = np.roll(target, y_roll, axis = -1)

        return torch.as_tensor(sample, dtype=torch.float), torch.as_tensor(target, dtype=torch.float)


if __name__ == '__main__':
    
    start = time.time()
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info(start)

    logging.info(torch.cuda.get_device_name(0))    
    num_cpus = os.cpu_count()
    logging.info("Number of CPUs available: {}".format(num_cpus)) 
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml_config", default='SFNO_exo.yaml', type=str)
    parser.add_argument("--config", default='sfno_backbone', type=str)
    args = parser.parse_args()
    
    params = YParams(os.path.abspath(args.yaml_config), args.config)

    world_size = int(os.environ.get('WORLD_SIZE', 1))
    if world_size > 1:
        dist.init_process_group(backend='nccl',
                            init_method='env://')
        local_rank = int(os.environ["LOCAL_RANK"])
        torch.cuda.set_device(local_rank)
        torch.backends.cudnn.benchmark = True
        rank = dist.get_rank()
    else:
        local_rank = 0
        rank = 0

    logging.info("World size: {}, local rank: {}".format(world_size, local_rank))

    path1 = "data/omega/"
    path2 = "data/omega_0_5/"
    path3 = "data/omega_4/"

    Y, X = np.meshgrid(np.linspace(-90, 90, 90), np.linspace(0, 361, 144))
    Yp, Xp = np.meshgrid(np.linspace(0, 4, 4), np.linspace(-90, 90, 90))

    dataset1 = GetDataset(path1, 1)
    dataset2 = GetDataset(path2, 5)
    dataset3 = GetDataset(path3, 10)

    sampler1 = DistributedSampler(dataset1) if world_size > 1 else None
    sampler2 = DistributedSampler(dataset2) if world_size > 1 else None
    sampler3 = DistributedSampler(dataset3) if world_size > 1 else None
    
    dataloader1 = DataLoader(dataset1,
            batch_size=params.batch_size,
            sampler=sampler1,
            shuffle=False,
            num_workers=2,
            pin_memory=torch.cuda.is_available(),
            drop_last=True
            )

    dataloader2 = DataLoader(dataset2,
            batch_size=params.batch_size,
            sampler=sampler2,
            shuffle=False,
            num_workers=2,
            pin_memory=torch.cuda.is_available(),
            drop_last=True
            )
   
    dataloader3 = DataLoader(dataset3,
            batch_size = params.batch_size,
            sampler=sampler3,
            shuffle=False,
            num_workers=2,
            pin_memory=torch.cuda.is_available(),
            drop_last=True
            )

    device = torch.cuda.current_device() if torch.cuda.is_available() else 'cpu'

    in_channels = params.in_channels
    out_channels = params.out_channels

    model = SFNO(params={},spectral_transform="sht",filter_type='linear',img_shape=(90, 144), in_chans=in_channels,
                                         out_chans=out_channels,
                                        num_layers=8,operator_type='dhconv',scale_factor=1,spectral_layers=3, embed_dim=256).to(device)
    
    if world_size > 1:
        model = DistributedDataParallel(model, broadcast_buffers=False )
    
    #logging.info("number of trainable parameters: ", count_parameters())
    optimizer = torch.optim.Adam(model.parameters(), lr = params.lr)
    loss_obj = LpLoss()
    criterion = nn.MSELoss()
    iters = 0
    tr_time = 0
    data_time = 0
    
    model.train()

    logging.info("length of dataloader1:{}".format(len(dataloader1)))
    logging.info("length of dataloader2:{}".format(len(dataloader2)))
    logging.info("length of dataloader3:{}".format(len(dataloader3)))

    logging.info(params.max_epochs)
    #restore_checkpoint('sfno_ckpt_amp_3.tar', model)
    
#     1 for 1, 5 for 0.5

    RR_1 = torch.full((params.batch_size,1,90,144), 1).to(device, dtype = torch.float)
    RR_2 = torch.full((params.batch_size,1,90,144), 5).to(device, dtype = torch.float)
    RR_3 = torch.full((params.batch_size,1,90,144),10).to(device, dtype = torch.float)

    scheduler = CosineAnnealingLR(optimizer, T_max=params.max_epochs)
   
    #torch.autograd.set_detect_anomaly(True)

    gscaler = torch.amp.GradScaler()

    for j in range(params.max_epochs):
        logging.info(j)
        i = 0
        if world_size > 1:
            sampler1.set_epoch(j)
            sampler2.set_epoch(j)
            sampler3.set_epoch(j)
        for (data1, data2, data3) in zip(dataloader1, dataloader2, dataloader3):

            inp1, tar1 = map(lambda x: x.to(device, dtype=torch.float), data1)
            inp2, tar2 = map(lambda x: x.to(device, dtype=torch.float), data2)
            inp3, tar3 = map(lambda x: x.to(device, dtype=torch.float), data3)

            inp1 = inp1.to(memory_format=torch.channels_last)
            tar1 = tar1.to(memory_format=torch.channels_last)
            inp2 = inp2.to(memory_format=torch.channels_last)
            tar2 = tar2.to(memory_format=torch.channels_last)
            inp3 = inp3.to(memory_format=torch.channels_last)
            tar3 = tar3.to(memory_format=torch.channels_last)
            
            model.zero_grad()

            with torch.amp.autocast(device_type='cuda', enabled=True):

                gen1 = model(inp1).to(device, dtype = torch.float)
                loss1 = loss_obj(gen1, tar1)

                gen2 = model(inp2).to(device, dtype = torch.float)
                loss2 = loss_obj(gen2, tar2)

                gen3 = model(inp3).to(device, dtype = torch.float)
                loss3 = loss_obj(gen3, tar3)

                if torch.isnan(gen1).any():# or torch.isnan(gen2).any():
                    logging.info("NaN detected in gen")
                    continue  # Skip the current iteration

            
            loss = loss1 + loss2 + loss3

            gscaler.scale(loss).backward()

            gscaler.step(optimizer)

            gscaler.update()

            #loss.backward()
            #optimizer.step()

            if torch.isnan(loss).any():
                logging.info("NaN detected in loss")
                continue  


            if i%10 == 0:
                if local_rank==0:
                    logging.info("save")
                    logging.info(f"Epoch {j}, iteration {i}, Loss 1: {criterion(gen1, tar1).item()}, Loss 2: {criterion(gen2, tar2).item()}, Loss 3: {criterion(gen3, tar3).item}") 
                    save_checkpoint(model)
                    logging.info(time.time() - start)

            i += 1

        scheduler.step()
            
    end = time.time()
    logging.info(end - start)

        
