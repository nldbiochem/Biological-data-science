import csv
import sys
from itertools import izip_longest, cycle
from collections import defaultdict
from collections import OrderedDict

def checkpoint(p, Queue_score):
    p_sum=sum(p)
    p_control_sum=sum(p[4:])
    if len(p)<5: #5 is where controls start
        if p_sum+(4-len(p))>=Queue_score: #4 equals diseased counts
            return 1
        else:
            return 0
    else:
        if p_sum-(2*p_control_sum)>=Queue_score:
            return 1
        else:
            return 0

def genes_with_pattern(up, down, data,diseased):
    count=sum(up)+sum(down)
    print(count)
    total=0
    locations_of_patterns=[]
    genes_locations=[]
    i=0
    for k in data:
        total=0
        if i<=diseased:
                  
          for j in range(0, len(up)):
            if up[j]==0 and down[j]==0 or data[k][j]==0:
                continue
                
            elif up[j]==data[k][j] or down[j]==-1*data[k][j]:
                
                total=1+total
                continue
                    
            else:
                continue
          if total==count:
              print('data is', data[k])
              locations_of_patterns.append('diseased'+k)
        i=i+1
        elif i>diseased:
            for j in range(0, len(up)):
                if up[j]==0 and down[j]==0 or data[k][j]==0:
                    continue
                
                elif up[j]==data[k][j] or down[j]==-1*data[k][j]:
                
                    total=1+total
                    continue
                    
            else:
                continue
        if total==count:
            locations_of_patterns.append('controls'+k)
    return locations_of_patterns

def calc_j(u, d, data, C, N,G, Queue_score):
    u_sum = sum(u)
    d_sum = sum(d)
    p = []
    i = 0
    for k,v in data.items():

        # calc a
        a = 0
        for i in range(G):
            if v[i]==1:
                a += v[i] * u[i]
        if a != u_sum:
            p.append(0)
            if checkpoint(p,Queue_score)==1:
                continue
            else:
                return 0, 0, 0


        # calc c
        c = 0
        for i in range(G):
            if v[i]==-1:
                c += v[i] * d[i] * -1
        if c != d_sum:
            p.append(0)

        else:
            p.append(1)

        if checkpoint(p, Queue_score) == 1:
            continue
        else:
            return 0, 0, 0
    # calc j
    diseased_p=0
    normal_p=0
    for i in range(0, C):
        diseased_p += float(p[i])
    float(diseased_p)
    diseased_J=diseased_p/C
    for i in range(C, C+N):
        normal_p+=float(p[i])
    normal_J=float( normal_p/N)
    J=float(diseased_J-normal_J)
    Queue_score=J
    
    return J, u, d


def file_processing(f):
   
    with open(f,'r') as csvf:
        dialect = csv.Sniffer().sniff(csvf.readline()) # finds the delimiters automatically
        csvf.seek(0)
                    # read file with dialect
        rdlistcsv = csv.reader(csvf,dialect)
                            # save to list of rows
        rowslist  = [list(filter(None,line)) for line in rdlistcsv]
        header = rowslist[0]
        data = {}
        for i,key in enumerate(header):
            ilist = [row[i] for row in rowslist if row[i] != key]
            data.update({key: ilist})

    for lists in data.values():
        for i in range(len(lists)):
            if  lists[i] == '-2':
                lists[i]=int(0)
            elif lists[i]== '-1':
                lists[i]=int(-1)
            elif lists[i] == '0':
                lists[i] = int(0)
            elif lists[i] == '1':
                lists[i]=int(1)
   
    for k,v in data.items():
        if type(v[0])==str:
            if type(v[0]) == str:
               del data[k]

                

    data_rev_sort=OrderedDict(sorted(reversed(data.items())))


    return data_rev_sort


def gene_nummbers(u,d):
    j=0
    up_list=[]
    d_list=[]
    for i in u:
        if i==1:
            up_list.append('High'+str(j))
        j= j+1
    j=0
    for i in d:
        if i ==1:
            d_list.append('Low'+str(j))
    
    return up_list, d_list

def main():
    df= sys.argv[1]
    Genes= sys.argv[2]
    Genes=int(Genes)
    print(Genes)

    processed_df=file_processing(df)

    D=sys.argv[3]#Diseased cases
    D=int(D)
    C=sys.argv[4]#number of controls
    output=sys.argv[5]
    C=int(C)
    Gstring=str(Genes)
    Queue_score=float('-inf')
    u=[]
    d1=[]
    d2=[]
    solution_u=[]
    solution_d=[]
    pattern_sum=0
    max_J=float('-inf')
    for i in range(0, 2**Genes):
        binarystr='{0:0'+Gstring+'b}'
        u_prime=list(binarystr.format(i))
        u=[int(i) for i in u_prime]
 
        for j in range(0, 2**Genes):
            d_prime=list(binarystr.format(j))
            d=[int(k) for k in d_prime]
            
            for l  in range(0, len(d)):
                if d[l] == 1 and u[l] == 1:
                    d[l]= 0    
           
            J,up,down=calc_j(u,d, processed_df,D,C,Genes,max_J)
            if J>=max_J and J>0:
                
                max_J=J
                if sum(up)+sum(down)>sum(solution_u)+sum(solution_d):
                    solution_u=up
                    solution_d=down
    pattern_locations=genes_with_pattern(solution_u,solution_d,processed_df,D)
    print('J is: ', max_J, ' up list is: ', solution_u, ' down list is: ',
            solution_d, pattern_locations)
    upl, dl= gene_nummbers(solution_u, solution_d)
    f=open(output, 'w')
    f.write('expression pattern'+'\n')
    f.write(str(upl)+'\n')
    f.write(str(dl)+'\n')
    f.write('cases and controls with pattern'+'\n')
    f.write(str(pattern_locations)+'\n')
    f.write('J score'+'\n')
    f.write(str(max_J))

    # binary processing
   # print(my_dict) ###/
if __name__ == "__main__":
    main()
