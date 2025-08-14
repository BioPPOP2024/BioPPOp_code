import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

Timez=pd.Series(range(0,11))
Timez=Timez*10

plt.plot(Timez,Timez)
plt.xlabel('Normalised Movement Time (%)', fontsize= 15)
plt.xticks(fontsize=12)
plt.show()