## Pacote cujo foco é ter ferramentas úteis para o estudo da distribuição de massa em galáxias.
## Aqui teremos desde formular para converter fluxo para massa, até formuladas com funções de GALFIT. 
## Tudo isso com o foco em estudar a distribuição de massa em galáxias do universo local. 
## Author: Yasmin Coelho
## Last update: 30/05/2022

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

### ATENÇÃO PROFESSOR, ESSAS FUNÇÕES ABAIXO NÃO SÃO USADAS NO TRABALHO DO MÓDULO 1 MAS EU JÁ TINHA ESCRITO ELAS NESSE DOCUMENTO E FIQUEI
### COM DÓ DE APAGAR ENTÃO IGNORA ELAS POR FAVOR. 
### ------ Essa função irá fazer a conversão de magnitude para fluxo. 
### magAB para fluxo em Jansky (Jy)

def mag_to_flux(mag):
	'''
	Sobre essa função: mag_to_flux(mag), onde:
	mag -- magnitude do objeto em mag AB
	Essa função irá fazer a conversão de magnitude para fluxo. magAB para fluxo em Jansky (Jy).
	'''

	flux = 10**((mag)*(-2)/5)*3631
	
	return flux
	
	
### ------ Essa função irá converter fluxo em massa
### Para isso, é necessario informar o fluxo do objeto em 3.6 e 4.5 micrometros e a distância 
### A fórmula de utilizada -- Eskew et al. 2012. 

def flux_to_mass(flux3,flux4,dist): ##Depois posso melhorar aqui para ser possível colcar a magnitude e ele converter aqui dentro também, antes de aplicar na fórmula
	'''
	Sobre essa função: flux_to_mass(flux3,flux4,dist), onde:
	PROMETO QUE NUM FUTURO PRÓXIMO EU EXPLICO. mas não hoje... perdão.
	'''
	mass = 10**(5.65)*(flux3)**(2.85)*(flux4)**(-1.85)*((dist)/(0.05))**2 
	
	return math.log10(mass)  #Devolve a massa em log10

###	AQUI COMEÇA AS FUNÇÕES QUE DE FATO UTILIZO NO MEU TRABALHO DO MÓDULO 1 (: 

### ------------
def intferrer(r,r0,s0,a=2.0,b=0.0):
    '''
    Ferrer Profile -- Essa função retorna o valor de fluxo da barra até determinado raio.
    Sobre essa função: intferrer(r,r0,s0,a=2.0,b=0.0), onde:
    r -- até onde será feita a integração para obter o fluxo da barra.  
    r0 -- raio de truncamento da barra [pix]
    s0 -- Fluxo da barra, convertido do brilho superficial. Unidade de fluxo/arcsec^2
    a e b são constantes nesse caso, valores adotados na decomposição das galáxias do S4G em Salo et al . 2015.
    
    Para mais informações sobre a função ferrer:
    Manual GALFIT -- https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf
    Salo et al. 2015 -- https://ui.adsabs.harvard.edu/abs/2015ApJS..219....4S/abstract
    '''
    return (2*pi*s0/(6*r0**4) * (r**6 - 3*r**4*r0**2 + 3*r**2*r0**4) )

### ------------
def b(n):
    '''
    Constante de normalização. Essa função existe para facilitar algumas contas. 
    '''
    return gammaincinv(2*n, 0.5)

### ------------
def exp_profile(D_s0, r, hr): 
    '''''
    Exponential Profile -- essa função retorna o fluxo do disco integrado até o infinito. 
    exp_profile(D_s0, r, hr), onde:
    D_s0 -- Fluxo do disco, convertido do brilho superficial. Unidade de fluxo/arcsec^2
    r -- até onde será feita a integração para obter o fluxo do disco.
    hr -- disk scale length [pix]   
    
    Para mais informações sobre a função exponential profile:
    Manual GALFIT -- https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf
    Salo et al. 2015 -- https://ui.adsabs.harvard.edu/abs/2015ApJS..219....4S/abstract
    '''''    
    S_r = D_s0*np.exp((-r)/hr) ##Galfit
    
    F_disk = 2 *pi *hr**2 * D_s0 
    
    return F_disk

### ------------
def exp_profile_enc(D_s0, r, hr, re):
    '''''
    Exponential Profile -- essa função retorna o fluxo do disco até certo raio.
    exp_profile_enc(D_s0, r, hr, re), onde:
    D_s0 -- Fluxo do disco, convertido do brilho superficial. Unidade de fluxo/arcsec^2
    r -- até onde será feita a integração para obter o fluxo do disco.
    hr -- disk scale length [pix]
    re -- raio efetivo do disco. re = hr  * 1.678.   [pix]
    
    Para mais informações sobre a função exponential profile:
    Manual GALFIT -- https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf
    Salo et al. 2015 -- https://ui.adsabs.harvard.edu/abs/2015ApJS..219....4S/abstract
    ''''' 
    D_n = 1.01  #Esse valor é 1, porém estava dando erro apenas com 1. 
    x = b(D_n) * (r/re)**(1.0/D_n)
    
    return exp_profile(D_s0, r, hr) * gammainc(2*D_n, x)
    
    
#### ------------------------------

def momento(data):
    '''
    Função que devolve o fotocentro da imagem utilizando o método do momento.
    Para que ela funcione da melhor forma possível é necessário que a imagem esteja cortada, englobando apenas o objeto de interesse.
    momento(data), onde:
    data -- valores dos pixels da imagem com o objeto.
    
    Essa função foi criada por Bruno Morgado.
    '''
    I_i = np.sum(data, axis=0)
    I_j = np.sum(data, axis=1)
    Ii_mean = np.sum(I_i)/len(I_i)
    Ij_mean = np.sum(I_j)/len(I_j)

    x_i = np.arange(data.shape[1])
    y_j = np.arange(data.shape[0])
    mask_i = (I_i - Ii_mean) > 0
    mask_j = (I_j - Ij_mean) > 0

    xc = np.sum((I_i - Ii_mean)[mask_i]*x_i[mask_i])/np.sum((I_i - Ii_mean)[mask_i])
    yc = np.sum((I_j - Ij_mean)[mask_j]*y_j[mask_j])/np.sum((I_j - Ij_mean)[mask_j])
    return xc, yc
    

#--------------------------------

def sigma_clip(value, sigma):
    '''
    Essa função analisa uma distribuição e retira delas os valores que estiverem fora do N sigma determinado na entrada. 
    A função permanece em loop até que só haja valores dentro do N sigma pedido.
    
    sigma_clip(value, sigma), onde:
    value -- valores da distribuição que será analisada
    sigma -- valor N sigma. ex: 1 sigma = 68.3% da distribuição; 2 sigmas = 95.4% da distribuição; 
                                3 sigmas = 99.7% da distribuição.
    '''
    sobra = value[np.absolute(value - value.mean()) > sigma*value.std(ddof=1)]
    while len(sobra) != 0:
        if len(sobra) != 0:
            s_clip = value[np.absolute(value - value.mean()) < sigma*value.std(ddof=1)]
            sobra = s_clip[np.absolute(s_clip - s_clip.mean()) > sigma*s_clip.std(ddof=1)]
            value = s_clip
            
    value_cut = value[np.absolute(value - value.mean()) < sigma*value.std(ddof=1)]
    sobra = value_cut[np.absolute(value_cut - value_cut.mean()) > sigma*value_cut.std(ddof=1)]
            
    return value_cut
    
    
#----------------------------------

