import sys
import os
import pandas as pd
import subprocess
import re 

rootdir = '/home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT' 
subdirs = [x for x in os.listdir(rootdir) if re.match(r'C',x)] 
#log_files = [x for x in os.listdir(rootdir) if re.match(r'pear_.*\.log',x)]

df = pd.DataFrame(columns=['sample','total_reads','uniquely_mapped','multiple_mapped','unmapped_reads'])
it = 0

#for each subdir with STAR alignment results, write commands to extract wanted numbers from Log.final.out file
#Take uniquely mapped and total number are in one line, multiple mapped and unmapped are in several lines
for x in subdirs:

    path_log_file=os.path.join(rootdir,x,x + '_Log.final.out')
    #build commands to give to subprocess
    command_total_reads="grep \"Number\sof\sinput\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9]*).*/$1/'"
    command_uniquely_mapped="grep \"Uniquely\smapped\sreads\snumber\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9]*).*/$1/'"
    command_multi_mapped_1="grep \"Number\sof\sreads\smapped\sto\smultiple\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9]*).*/$1/'  "
    command_multi_mapped_2="grep \"Number\sof\sreads\smapped\sto\stoo\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9]*).*/$1/'"
    command_unmapped_1="grep \"\%\sof\sreads\sunmapped\:\stoo\sshort\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9.]*).*/$1/'"
    command_unmapped_2="grep \"\%\sof\sreads\sunmapped\:\stoo\smany\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9.]*).*/$1/'"
    command_unmapped_3="grep \"\%\sof\sreads\sunmapped\:\sother\" " + path_log_file + " | perl -pe 's/.*\|\s*([0-9.]*).*/$1/'"
 
   #run subprocess commands 
    total_reads = int(subprocess.check_output(command_total_reads , shell=True).strip('\n'))
    uniquely_mapped = int(subprocess.check_output(command_uniquely_mapped , shell=True).strip('\n'))
    multi_mapped_1 = int(subprocess.check_output(command_multi_mapped_1 , shell=True).strip('\n'))
    multi_mapped_2 = int(subprocess.check_output(command_multi_mapped_2 , shell=True).strip('\n'))
    unmapped_1 = float(subprocess.check_output(command_unmapped_1 , shell=True).strip('\n'))
    unmapped_2 = float(subprocess.check_output(command_unmapped_2 , shell=True).strip('\n'))
    unmapped_3 = float(subprocess.check_output(command_unmapped_3 , shell=True).strip('\n'))
        
    #if numbers are all extracted , calculate unmapped and multiple mapped and write to dataframe 
    if total_reads and uniquely_mapped:

        multi_mapped = multi_mapped_1 + multi_mapped_2
        unmapped = int(total_reads*((unmapped_1+unmapped_2+unmapped_3)*0.01))
            
        new_row = [x,total_reads,uniquely_mapped,multi_mapped,unmapped]         #x is the subdir and identical with sample_id
        df.loc[it] = new_row 
        
        it += 1
        
        print x + ': worked'
    else:
        print x + ': something is wrong'

#save dataframe
df.to_csv('/home/stroemic/hiwi_16/analysis/STAR_alignment/5Seq_processed_IVT/5Seq_processed_IVT_STAR_statistics.csv', index = False)
