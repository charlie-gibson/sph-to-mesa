def sph_input_reader():
    # Create an empty dictionary to store the data
    sph_data = {}
    
    # Read in file
    with open('sph.input', 'r') as f:
        n = 0
        
        for line in f:
            try:
                parts = line.strip().split('!')
                value_part = parts[0].split('=')

                if len(value_part) >= 2:
                    variable_name = value_part[0].strip().lower()
                    if variable_name[0]=='n' or variable_name[0]=='N':
                        value=int(value_part[1].strip().replace(',',''))
                    else:
                        try:
                            value = float(value_part[1].strip().replace(',',''))
                        except:
                            value=value_part[1].strip().replace(',','').strip("'").strip('"')
                    
                    #descriptor = None
                    #if len(parts) == 2:
                    #    descriptor = parts[1].strip()
                    
                    #sph_data[variable_name] = {'value': value, 'descriptor': descriptor}
                    sph_data[variable_name] = value
                    n += 1
                else:
                    print(f"Ignoring invalid line {n}: {line}")
            except Exception as e:
                print(f"Error processing line {n}: {e}")

        for keys,vals in zip(sph_data.keys(),sph_data.values()):
            print(f'{keys:15} {vals}')
    
    return sph_data

