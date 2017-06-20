import sys


def interval_avg(depth_interval, interval_avg, sample_mean):
    with open(depth_interval, 'r') as interval_in:
        with open(interval_avg, 'w+') as interval_out:
            #skip first line in depth file
            interval_in.readline()
            
            for line in interval_in:
                #no need to trip carriage return as only extracting first and third column
                values = line.split('\t')
                
                interval = values[0]
                interval_cvg = values[2]
                normalize_value = float(interval_cvg)/float(sample_mean)
                interval_out.write(interval + '\t' + str(normalize_value) + '\n')
                    
            

def sample_mean_cvg(depth_summary):
    with open (depth_summary, 'r') as depth_file:
        #skip first line
        depth_file.readline()
        #get mean cvg from 2nd line
        depth_values = depth_file.readline()
        values = depth_values.split('\t')
        mean_cvg = values[3]
        return mean_cvg
        
        
        
        

def main():
    depth_interval = sys.argv[1]
    depth_summary = sys.argv[2]
    normalized_cvg = sys.argv[3]
    
    sample_mean = sample_mean_cvg(depth_summary)
    
    interval_avg(depth_interval, normalized_cvg, sample_mean)
main()    
