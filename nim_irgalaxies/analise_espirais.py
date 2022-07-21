## Pacote cujo foco é ter ferramentas úteis para o estudo da distribuição de massa em galáxias.
## Aqui teremos desde formular para converter fluxo para massa, até formuladas com funções de GALFIT. 
## Tudo isso com o foco em estudar a distribuição de massa em galáxias do universo local. 
## Author: Yasmin Coelho
## Last update: 21/07/2022

import numpy as np
import math
import numpy as np
from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D, Sersic2D
import scipy.integrate as integrate
from scipy.integrate import quad, dblquad
from numpy import pi, exp
from scipy.special import gamma, gammaincinv, gammainc
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u


#----------------------------------

def analise_espirais(column_type, plot=True):
    '''
    Essa função analisa uma lista de "Hubble morphological type" e verifica quantas galáxias são espirais(S), quantas 
    possuem uma barra(SB) e quantas barras são fortes (SAB).
    
    Essa função foi feita para ler uma coluna de um data frame, sendo necessário indicar a tabela e a coluna que lista 
    os tipos morfológicos.
    
    ----
    
    Temos  analise_espirais(column_type, plot=True), onde:
    column_type -- lista com os tipos morfológicos. tipo string. 
    plot = True -- devolve o histograma de tipo morfolófico para a lista de entrada
    
    ex: analise_espirais(column_type = s4g_2267.type, plot=True)
    '''
    
    lista_type = np.asarray(column_type, dtype=str)

    lista_s = [] #lista para espirais
    lista_sb = [] #lista para espirais barradas
    lista_sab = [] #lista para espirais com barras fortes
    for i in range(len(lista_type)):
    
        if 'S' in lista_type[i]:
            lista_s.append(lista_type[i])
      
        if 'SB' in lista_type[i] or 'SAB' in lista_type[i]:
            lista_sb.append(lista_type[i])
        
        if 'SAB' in lista_type[i]:
            lista_sab.append(lista_type[i])

    print('Amostra de', len(lista_type) ,'galáxias\n')
    print('Em nossa amostra de: ', len(lista_type), 'galáxias temos: ', len(lista_s), 'galáxias espirais')
    print('onde: ', len(lista_sb), 'possuem uma barra')
    print('sendo : ', len(lista_sab), 'são barras fortes')
    print('\nPorcentagem de galaxias espirais barradas:  %.3f' %(len(lista_sb)*100/len(lista_type)))
    print('Porcentagem de galaxias espirais com barras fortes:  %.3f' %(len(lista_sab)*100/len(lista_type)))
    print('\n')
    if plot:
        fig = plt.figure(figsize=(16,5))
        ax1 = fig.add_subplot(1,1,1)

        ax1.hist(sorted(lista_type),bins=30)
        ax1.set_title('Analise de tipo morfológico de uma amostra com tamanho: %.f' %float(len(lista_type))) 
   
        plt.show()
