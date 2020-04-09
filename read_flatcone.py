import numpy as np

# reads in flatcone data as a dictionary
# no other reduction of the data is performed
def read_flatcone(filename):
    with open(filename) as f:
        raw = f.readlines()

    # get indices for start of data/multi
    for i, line in enumerate(raw):
        if 'DATA_:' in line:
            data_start = i
        elif 'MULTI:' in line:
            multi_start = i

    # get parameters needed
    params = ['AS', 'BS', 'CS', 'AA', 'BB', 'CC', 'AX', 'AY', 'AZ', 'CHAN', 'KFIX', 'DM', 'DA', 'SM', 'SA', 'SS', 'FX']
    param_dict = {}
    for line in raw[:data_start]:
        ls = line.replace(',','').split()
        for param in params:
            p = param + '='
            if p in ls:
                index = ls.index(p) + 1
                param_dict[param] = float(ls[index])
            
    # extract data and headers
    header = raw[data_start+1].split()
    data_str = raw[data_start+2:multi_start]
    data_str = [x.split() for x in data_str]

    # convert each line of data to floats
    data = []
    for d in data_str:
        data.append([float(x) for x in d])

    # conver to numpy array
    data = np.array(data)

    # extract flatcone data
    multi_str = raw[multi_start+1:]
    multi_str = [x.split() for x in multi_str]

    # convert each line of multi to floats
    multi = []
    for m in multi_str:
        multi.append([float(x) for x in m])

    # convert to numpy array
    multi = np.array(multi)

    # save as dictionary
    data_dict = {}
    data_dict['multi'] = multi
    data_dict['params'] = param_dict

    for i, h in enumerate(header):
        data_dict[h] = data[:,i]

    return data_dict

