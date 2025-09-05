# from _typeshed import NoneType
import pdb
import torch
from torch import nn
# from matplotlib.pyplot import imread
import numpy as np
import os
import torchvision.transforms as transforms
import torch.nn.functional as F


from efficientnet_pytorch import EfficientNet

mse_loss = torch.nn.MSELoss(reduction='none')

def loss_function(y, y_pred):
    return mse_loss(y, y_pred).sum(axis=1)


class FocalLoss(nn.Module): 
    def __init__(self, alpha=torch.tensor([0.22484802, 0.10919441, 0.05926421, 0.18060437, 0.19005628,
       0.18120403, 0.05482868]), gamma=2.0, reduction='mean'):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction

    def forward(self, inputs, targets):
        BCE_loss = F.cross_entropy(inputs, targets, reduction='none')
        pt = torch.exp(-BCE_loss)
        if self.alpha is not None:
            alpha = self.alpha[targets]
        else:
            alpha = 1
        F_loss = alpha * (1 - pt) ** self.gamma * BCE_loss

        if self.reduction == 'mean':
            return torch.mean(F_loss)
        elif self.reduction == 'sum':
            return torch.sum(F_loss)
        else:
            return F_loss


class MLP2(nn.Module):
    '''
    Multilayer Perceptron.
    '''
    def __init__(self):
        super(MLP, self).__init__()

        # Layers
        self.linear1 = nn.Linear(6, 32)
        self.linear2 = nn.Linear(32, 64)
        self.linear3 = nn.Linear(64, 32)
        self.linear4 = nn.Linear(32, 7)

        # ReLU Activation
        self.relu = nn.ReLU()

        # Loss function
        self.criterion_gaze = FocalLoss()

        # Initialize weights and biases
        nn.init.xavier_uniform_(self.linear1.weight)
        nn.init.zeros_(self.linear1.bias)
        nn.init.xavier_uniform_(self.linear2.weight)
        nn.init.zeros_(self.linear2.bias)
        nn.init.xavier_uniform_(self.linear3.weight)
        nn.init.zeros_(self.linear3.bias)
        nn.init.xavier_uniform_(self.linear4.weight)
        nn.init.zeros_(self.linear4.bias)

    def forward(self, x):
        '''Forward pass'''
        x = self.relu(self.linear1(x))
        x = self.relu(self.linear2(x))
        x = self.relu(self.linear3(x))
        x = self.linear4(x)
        return x

    def predict(self, x):
        '''Prediction step, returns class indices'''
        x_pred = self.forward(x)
        return torch.argmax(x_pred, dim=1)

    def loss(self, x, y):
        '''Computes loss between predicted and true labels'''
        y_pred = self.forward(x)
        return self.criterion_gaze(y_pred, y)
  
class DataAugmentationLayer(nn.Module):
    def __init__(self, mean, std, noise_std_factor=0.05):
        """
        Layer per applicare data augmentation al volo.

        Args:
            mean: torch.Tensor, media delle feature.
            std: torch.Tensor, deviazione standard delle feature.
            noise_std_factor: float, fattore di scala per il rumore.
        """
        super(DataAugmentationLayer, self).__init__()
        self.mean = mean
        self.std = std
        self.noise_std_factor = noise_std_factor  # Controlla l'entità del rumore

    def forward(self, x):
        if self.training:  # Applica la data augmentation solo in training
            #print("Before augmentation: ", {x[0]})
            noise = torch.randn_like(x) * (self.std * self.noise_std_factor)
            x = x + noise
            #print("After augmentation: ", {x[0]})
        return x


class MLP(nn.Module):
    '''
    Multilayer Perceptron.
    '''
    def __init__(self, mean, std):
        super().__init__()
        #self.criterion_gaze = nn.CrossEntropyLoss()
        self.augment = DataAugmentationLayer(mean, std, noise_std_factor=0.3)
        self.criterion_gaze = FocalLoss()
        self.layers = nn.Sequential(
          nn.Linear(6, 32),
          nn.LeakyReLU(0.1),
          nn.Dropout(0.3),
    
          nn.Linear(32, 64),
          nn.LeakyReLU(0.1),
          nn.Dropout(0.3),
    
          nn.Linear(64, 32),
          nn.LeakyReLU(0.1),
          nn.Dropout(0.3),
    
          nn.Linear(32, 7),
          nn.Softmax(dim=1)  # Softmax per garantire la distribuzione delle probabilità
        )


    '''def forward(self, x):
        Forward pass
        return self.layers(x)'''

    def forward(self, x):
        activations = []
        x = self.augment(x)
        for layer in self.layers:
            x = layer(x)
            if isinstance(layer, (nn.LeakyReLU, nn.Softmax)):
                activations.append(x.clone().detach())
        return x, activations

    def predict(self, x):
        output, _ = self.forward(x)
        return torch.argmax(output, dim = 1)

    def loss(self, x, y):
        y_pred, _ = self.forward(x)
        return self.criterion_gaze(y_pred, y)
'''
    def predict(self, x):
        
        x_pred = self.layers(x)
        # x_final = torch.zeros(7, device=x_pred.device)
        # x_final[torch.argmax(x_pred)] = 1
        # return x_final
        return torch.argmax(x_pred, dim=1)

    def loss (self, x,y):
        # y_pred = self.predict(x)
        y_pred = self.forward(x)


        loss = self.criterion_gaze(y_pred, y)
        return loss'''

######################################
class CNN(nn.Module):
    '''
    Convolutional Neural Network - based on Efficient Net and WheNet
    '''
    def __init__(self):
        super().__init__()
        self.criterion_immagini = nn.CrossEntropyLoss()
        self.criterion_gaze = nn.MSELoss()
        #self.criterion_gaze = nn.CrossEntropyLoss()
        self.efficientnet =  torch.hub.load('NVIDIA/DeepLearningExamples:torchhub', 'nvidia_efficientnet_b0', pretrained=True)
        self.efficientnet.eval().to(torch.device("cuda"))
        # utils = torch.hub.load('NVIDIA/DeepLearningExamples:torchhub', 'nvidia_convnets_processing_utils')
        self.avgpool = nn.Sequential(
            nn.Flatten(),
            nn.Linear(1000, 500),
            nn.Softmax(dim=1),
            nn.Linear(500, 32),
            nn.ReLU(),
            nn.Linear(32, 6)
        )
        # self.denselay = nn.Linear(15, 6)
        # self.softmax = nn.Softmax()

    def forward(self, x):
        '''Forward pass'''

        # newx = torch.zeros(x.shape[0], x.shape[3], x.shape[1], x.shape[2]).cuda()
        # newx[:,0,:,:] = x[:,:,:,0]
        # newx[:,1,:,:] = x[:,:,:,1]
        # newx[:,2,:,:] = x[:,:,:,2]
        
        # newx_normalize= tfms(newx)
        # 
        prex = self.efficientnet(x)
        #prexnew = prex
        prexnew = prex.detach()
        #import pdb; pdb.set_trace()
        y = self.avgpool(prexnew)
        
        return y
    def predict(self, x):
        '''Forward pass'''
        x_pred = self.forward(x)
        # x_prob = self.softmax(x_pred)
        x_final = torch.zeros(6)
        x_final[torch.argmax(x_pred)] = 1
        return x_final

    def loss (self, x,y):
        y_pred = self.forward(x)
        # import pdb; pdb.set_trace()
        #ytrue =((y == 1).nonzero(as_tuple=True)[1])
        ytrue = y
        loss = self.criterion_gaze(y_pred, ytrue)
        return loss

class MLP_CNN(nn.Module):
    '''
    Multilayer Perceptron and Convolutional Neural Network
    '''
    def __init__(self, model_mlp = 'None', model_cnn ='None'):
        super().__init__()
        self.criterion_all = nn.CrossEntropyLoss()
        self.efficientnet =  model_cnn
        self.mlp =  model_mlp
        # utils = torch.hub.load('NVIDIA/DeepLearningExamples:torchhub', 'nvidia_convnets_processing_utils')
        self.avgpool = nn.Linear(12,6)
        # self.denselay = nn.Linear(15, 6)
        # self.softmax = nn.Softmax()

    def forward(self,x):
        '''Forward pass'''

        gaze, immagine = x
        feat_immagine = self.efficientnet(immagine)
        feat_gaze = self.mlp(gaze)
        total_x =  torch.cat([feat_immagine, feat_gaze],1)
        y = self.avgpool(total_x)
        
        return y
    def predict(self, x):
        '''Forward pass'''
        x_pred = self.forward(x)
        # x_prob = self.softmax(x_pred)
        x_final = torch.zeros(6)
        x_final[torch.argmax(x_pred)] = 1
        return x_final

    def loss (self, x,y):
        y_pred = self.forward(x)
        ytrue =((y == 1).nonzero(as_tuple=True)[1])
        
        loss = self.criterion_all(y_pred, ytrue)
        return loss