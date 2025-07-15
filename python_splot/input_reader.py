def sph_input_reader():
    # Create an empty dictionary to store the data
    sph_data = {}
    
    # Read in sph.input
    with open('sph.input', 'r') as f:
        n = 0
        
        for line in f:
            try:
                parts = line.strip().split('!')
                value_part = parts[0].split('=')

                if len(value_part) >= 2:
                    variable_name = value_part[0].strip().lower()
                    # all variables beginning with n are integers
                    if variable_name[0]=='n' or variable_name[0]=='N':
                        value=int(value_part[1].strip().replace(',',''))
                    else:
                        try:
                            # removes , from sph.input data
                            value = \
                                float(value_part[1].strip().replace(',',''))
                        except:
                            # ensures strings don't have two sets of
                            # quotes
                            value= \
                                value_part[1].strip().replace(',','')\
                                    .strip("'").strip('"')

                    # creates a key for the sph.input value
                    # and sets that value in the dictionary
                    # to the value assigned in sph.input
                    sph_data[variable_name] = value
                    n += 1
                else:
                    print(f"Ignoring invalid line {n}: {line}")
            except Exception as e:
                print(f"Error processing line {n}: {e}")

        # prints sph.input data
        for keys,vals in zip(sph_data.keys(), sph_data.values()):
            print(f'{keys:15} {vals}')

    return sph_data

