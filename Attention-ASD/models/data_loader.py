import torch
from torch import nn

import json

import os
import torch
from torch.utils.data import Dataset
import numpy as np
import pickle
from scipy import io as sp_io

import numpy as np
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import pickle
from copy import deepcopy
from scipy import io as sp_io
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.compose import ColumnTransformer
import torchvision.transforms as transforms
from PIL import Image
image_size = 224
tfms = transforms.Compose([transforms.Resize(image_size), transforms.CenterCrop(image_size), 
                           transforms.ToTensor(),
                           transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),])

class AttentionDataset(Dataset):
    # TODO: Missing data corruption
    def __init__(self, path, part='train', cut=np.inf, split = 'random'):
        self.part = part
        self.data_x = None
        self.data_y = None
        if part == 'train':
            # self.create_train_test_split(split)
            self.read_data('train')
        elif part == 'test':
            self.read_data('test')
        elif part == 'val':
            self.read_data('val')
    def _convert_one_hot_entry(self, array, idx):
        first_part = array[:idx]
        action = int(array[idx])
        second_part = array[idx+1:]

        one_hot_vector = np.zeros(4)
        one_hot_vector[action-1] = 1

        transformed_vector = np.concatenate([first_part, one_hot_vector, second_part], axis=0)

        return transformed_vector

    def read_data(self, label):

        # import pdb; pdb.set_trace()
        if label == 'train':
            trainAttrX, trainY1 = pickle.load(open('/content/drive/MyDrive/PW_2024/pickles/no_hbb_me_train.pkl', 'rb'), encoding='latin1')

            self.gaze = trainAttrX
            self.output = trainY1
            # import pdb; pdb.set_trace()
        if label == 'test':
            testAttrX, testY1 = pickle.load(open('/content/drive/MyDrive/PW_2024/pickles/no_hbb_me_test.pkl', 'rb'), encoding='latin1')

            self.gaze = testAttrX
            self.output = testY1
        if label == 'val':
            val_data, val_output_finale = pickle.load(open('/content/drive/MyDrive/PW_2024/pickles/no_hbb_me_val.pkl', 'rb'), encoding='latin1')

            self.gaze = val_data
            self.output = val_output_finale
            # import pdb; pdb.set_trace()

    def __len__(self):
        return len(self.gaze)

    def __getitem__(self, idx): 
        gaze = self.gaze[idx]
        output = self.output[idx]
        immagine = self.immagine[idx]
        return torch.tensor(gaze).float(), torch.tensor(immagine), torch.tensor(output).float()