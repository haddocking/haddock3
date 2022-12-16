import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import pdb



class cnn_class(nn.Module):

    def __init__(self,input_shape):
        # input_shape: (C, W, H, D)
        super(cnn_class,self).__init__()

        self.bn0 = nn.BatchNorm3d(input_shape[0])
        self.conv1 = nn.Conv3d(in_channels = input_shape[0], out_channels = 6,kernel_size=3)
        self.bn1 = nn.BatchNorm3d(6)
        self.mp1 = nn.MaxPool3d((3,3,3))
        self.conv2 = nn.Conv3d(in_channels = 6, out_channels = 6,kernel_size=3)   
        self.bn2 = nn.BatchNorm3d(6)
        self.mp2 = nn.MaxPool3d((3,3,3))
        size = self._get_conv_outputSize(input_shape)
        self.fc_1 = nn.Linear(in_features=size,out_features = 6)
        self.bn3 = nn.BatchNorm1d(6)
        self.fc_2 = nn.Linear(6,2)


    def _get_conv_outputSize(self,shape):
        num_data_points = 10
        inp = Variable(torch.rand(num_data_points,*shape))
        out = self._forward_features(inp)
        return out.data.view(num_data_points,-1).size(1)

    def _forward_features(self,x):
        x = F.max_pool3d(F.relu(self.bn1(self.conv1(x))),3)       
        x = F.max_pool3d(F.relu(self.bn2(self.conv2(x))),3)
        
        ''''

        x = F.relu(self.conv3(x))
        x = self.bn3(x)

        
        x = F.relu(self.conv4(x))
        x = self.pool4(x)
        x = self.bn4(x)

        x = F.relu(self.conv5(x))
        x = self.pool5(x)
        x = self.bn5(x)
        '''
        return x

    def forward(self,x):
        x = self._forward_features(x)
        x = x.view(x.size(0),-1)
        x = F.relu(self.fc_1(x))
        x = self.bn3(x)
        x = self.fc_2(x)
        return x
