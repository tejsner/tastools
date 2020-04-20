import numpy as np

# reads in flatcone data as a dictionary
# no other reduction of the data is performed
def read_flatcone(filename):
    with open(filename) as f:
        raw = f.readlines()

    # initialize dictionary to return
    r = {}
    r['INFO'] = {}
    multi_start = None # parsing still works for non-flatcone files

    # get indices for start of data/multi/parameters
    for i, line in enumerate(raw):
        if 'DATA_:' in line:
            data_start = i
        
        if 'MULTI:' in line:
            multi_start = i
        
        if line.startswith('VVVVV'):
            params_start = i

    # save all the parameters
    for line in raw[params_start+1:data_start]:
        # dictionary key is the 5 letter prefix in the header
        prefix = line[:5]

        # if the string contains '=', record the various parameters
        # otherwise save it as a string in r['info'][prefix]
        if '=' in line[7:]:            
            if prefix not in r:
                r[prefix] = {}
            
            # remove whitespace and newline, split by ','
            param_list = line[7:].replace(' ', '').replace('\n','').split(',')
            for param in param_list:
                # get key, value pairs by splitting on '='
                key, value = param.split('=')
                
                # convert to float if possible and save
                try:
                    r[prefix][key] = float(value)
                except:
                    r[prefix][key] = value
        else:
            r['INFO'][prefix] = line[7:].replace('\n', '')
            
    # extract data and headers
    header = raw[data_start+1].split()
    data_str = raw[data_start+2:multi_start]
    data_str = [x.split() for x in data_str]

    # convert each line of data to floats
    data = []
    for d in data_str:
        data.append([float(x) for x in d])

    # convert to numpy array
    data = np.array(data)

    # save each column
    r['DATA_'] = {}
    for i, h in enumerate(header):
        r['DATA_'][h] = data[:,i]

    if r['INFO']['TYPE_'] == 'flatcone':
        # extract flatcone data
        multi_str = raw[multi_start+1:]
        multi_str = [x.split() for x in multi_str]

        # convert each line of multi to floats
        multi = []
        for m in multi_str:
            multi.append([float(x) for x in m])

        # convert to numpy array and save
        r['MULTI'] = np.array(multi)

    return r
