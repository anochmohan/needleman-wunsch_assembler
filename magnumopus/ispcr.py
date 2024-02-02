#!/usr/bin/env python3

import subprocess
import tempfile

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:

    #primer_file and assembly_file are in q1.py

    cmd = ["blastn", "-task", "blastn-short", "-query", primer_file, "-subject" ,assembly_file, "-outfmt", "6 std qlen"] #blast command
    
    # subproccess runs the cmd command, and pipes the std.out
    result = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    sub_output = subprocess.Popen(('awk', '$3>80 && $4==$13'), stdin=result.stdout, stdout=subprocess.PIPE)
    output = subprocess.check_output(('sort', '-k 9', '-g'), stdin=sub_output.stdout, text=True) # one big string
    
    line_list = output.splitlines() # splits the output string at \n
    
    sorted_hits = []
    # This is done to get the list[list[str]] style
    for line in line_list:
        x = line.split("\t")
        sorted_hits.append(x)
        
    #return output


    # sorted hits and max_amplicon_size are in q2.py

    forward_primer = []
    reverse_primer = []
    
    #trying to separate the list into two, as forward and reverse
    for line in sorted_hits:
        if line[8] > line[9]:
            reverse_primer.append(line)
        else:
            forward_primer.append(line)
            
    hit_pairs = [] # only add primers that meet the criteria       
    for i in forward_primer:
        for j in reverse_primer:
            if int(j[8]) - int(i[9]) > max_amplicon_size:
                pass
            elif int(j[8]) - int(i[9]) <= 0:
                pass
            else:
                hit_pairs.append((i,j)) #stroe it as a tuple in a list
                
    #return hit_pairs

    bed_list = [] 
    # Trying to create a bed file style
    for a,b in hit_pairs:
        bed_list.append(a[1]+ "\t" + a[9] + "\t" + str(int(b[9])-1)) #a[1] = name, a[9] = start, b[9]-1 = end (-1 cuz seqtk is 0 indexed. blast is 1 indexed)
    bed_string = "\n".join(bed_list)
    
    #tempfile is used cuz seqtk input takes in a file, not string
    with tempfile.NamedTemporaryFile(mode="w+") as tmp:
        tmp.write(bed_string)
        tmp.seek(0) #sets cursor back at the start
        temp_filepath = tmp.name #stores temp filepath
        
        cmd = ["seqtk", "subseq", assembly_file, temp_filepath] #found this on seqtk giithub
        
        result = subprocess.run(cmd, capture_output=True, text=True)
    
    
    return result.stdout
