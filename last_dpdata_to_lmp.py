import dpdata as dp

data=dp.LabeledSystem('dpdata','deepmd/npy')
data_last=data.sub_system(-1)
data_last.to('lmp','conf.data_NO_MASS')


