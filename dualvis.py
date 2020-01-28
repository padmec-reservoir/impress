for i in range(M.vol_indicator.shape[1]):
    M.press[M.vol_indicator[:, i]]=i
M.pressure[:]=1
for i in range(M.pre_vol_indicator.shape[0]):
   if M.pre_vol_indicator[i].sum()==1:
       M.pressure[i]=2000
M.save_variables('visu')
