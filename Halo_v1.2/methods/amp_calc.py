import sys



def sample_mean_cvg(depth_summary):
    with open (depth_summary, 'r') as depth_file:
        #skip first line
        depth_file.readline()
        #get mean cvg from 2nd line
        depth_values = depth_file.readline()
        values = depth_values.split('\t')
        mean_cvg = values[2]
        return mean_cvg
        
#this method normalizes the sample to its own mean coverage     
def interval_normalized(depth_interval,sample_mean):
    with open(depth_interval, 'r') as interval_in:
            
            interval_in.readline()

            normalized_dict= {}
            #this method makes the assumption that there are the same number of lines
            # in your hospots normilization as there are in your intervals file
            for line in interval_in:
                #no need to trip carriage return as only extracting first and third column
                values_sample = line.split('\t')
                interval = values_sample[0]
                interval_cvg = values_sample[2]
                mean_normed = float(interval_cvg)/float(sample_mean)
                normalized_dict[interval] = mean_normed

            return normalized_dict
#this method finds the same intervals in a normals panel, and normalizes to this value
def normalize_to_normals(normals_intervals, sample_normed_dict):
    gene_dict = {}
    with open(normals_intervals, 'r') as normals_intervals:
        for line in normals_intervals:
            line=line.strip('\n')
            values = line.split('\t')
            norm_interval = values[0]
            norm_mean = values[1]
            norm_gene = values[2]
            norm_gene = norm_gene.strip("\r\n")
            if norm_interval in sample_normed_dict:
                sample_normed_normal = (sample_normed_dict.get(norm_interval) / float(norm_mean))
                if norm_gene not in gene_dict.keys():
                    gene_dict[norm_gene] = [sample_normed_normal]
                else:
                    gene_dict[norm_gene].append(sample_normed_normal)
    return gene_dict
                
            

def calc_write_avg(dictionary, out_file, sample_name):
    with open(out_file, 'w') as out:
        out.write('#SampleName\tGene\tAmp\n')
        for key,value in dictionary.iteritems():
            avg_value = sum(value)/len(value)
            out.write(sample_name + '\t' + key + '\t'  + str("%.2f" % avg_value) + '\n')





def main():
    depth_interval = sys.argv[1]
    depth_summary = sys.argv[2]
    normalized_cvg = sys.argv[3]
    out_file = sys.argv[4]
    sample_name = sys.argv[5]

    sample_mean = sample_mean_cvg(depth_summary)
    sample_dict = interval_normalized(depth_interval, sample_mean) 
    normalized_dict = normalize_to_normals(normalized_cvg, sample_dict) 
    print str(normalized_dict)
    calc_write_avg(normalized_dict, out_file, sample_name)
main()

