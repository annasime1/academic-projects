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

[]
class FocalLoss(nn.Module): 
    def __init__(self, alpha=torch.tensor([0.138, 0.033, 0.039, 0.620, 0.663, 2.310, 0.028]), gamma=2.0, reduction='mean'):
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


class MLP_FT(nn.Module):
    '''
    Multilayer Perceptron.
    '''
    def __init__(self, mean, std):
        super().__init__()
        # self.criterion_gaze = nn.CrossEntropyLoss()
        self.augment = DataAugmentationLayer(mean, std, noise_std_factor=0.3)
        self.criterion_gaze = FocalLoss()
        self.layers = nn.Sequential(
            nn.Linear(6, 32),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(32, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(32, 7), 
            nn.Softmax(dim=1) #aggiunta da noi per garantire la distribuzione come delle probabilità
        )

    '''def forward(self, x):
        Forward pass
        return self.layers(x)'''

    def forward(self, x):
        activations = []
        x = self.augment(x)
        for layer in self.layers:
            x = layer(x)
            if isinstance(layer, (nn.ReLU, nn.Softmax)):
                activations.append(x.clone().detach())
        return x, activations

    def predict(self, x):
        output, _ = self.forward(x)
        max_class = torch.max(output).item()
        return torch.argmax(output, dim=1)
        #return torch.argmax(output, dim = 1)

    def loss(self, x, y):
        y_pred, _ = self.forward(x)
        return self.criterion_gaze(y_pred, y)


