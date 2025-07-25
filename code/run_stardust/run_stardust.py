import os
import subprocess
import shutil
#from .build import build_start
from datetime import datetime
start_time = datetime.now()

def run():
        print("Setting up directory structure")
        if os.path.exists("Stardust_results"):
            shutil.rmtree("Stardust_results")
        os.mkdir("Stardust_results")
        os.mkdir("Stardust_results/build_output")
        os.mkdir("Stardust_results/visualization_output")
        for i in range(1,5):
            os.mkdir("Stardust_results/build_output/"+str(i)+"_pass")

        for i in range(1,5):
            os.mkdir("Stardust_results/visualization_output/"+str(i)+"_pass")

        dataset_path = input("Enter the dataset path\n")
        type = input("Enter the expression file type\n")

        for i in range(1,5):
            print("==============================================")
            print("PASS STARTED: ",i)
            print("==============================================")
            #build_start(dataset_path,"Stardust_results/build_output/"+str(i)+"_pass",type,i)
            if i==1:
                subprocess.call(['python','-W','ignore','stardust/run_stardust/build.py','-i',dataset_path,'-t',type,'-o','Stardust_results/build_output/'+str(i)+'_pass/','-n_pass',str(i)])
            else:
               subprocess.call(['python','-W','ignore','stardust/run_stardust/build.py','-i',dataset_path,'-t',type,'-o','Stardust_results/build_output/'+str(i)+'_pass/','-n_pass',str(i),'-pca_n','0'])
            if i==4:
                subprocess.call(['python','-W','ignore','stardust/run_stardust/embedding.py','-inp_data',dataset_path,'-i','Stardust_results/build_output/'+str(i)+'_pass/',
                             '-o','Stardust_results/visualization_output/'+str(i)+'_pass/','-n_pass',str(i)])
            else:
                subprocess.call(['python','-W','ignore','stardust/run_stardust/visualization.py','-inp_data',dataset_path,'-i','Stardust_results/build_output/'+str(i)+'_pass/',
                             '-o','Stardust_results/visualization_output/'+str(i)+'_pass/','-n_pass',str(i)])

        end_time = datetime.now()
        print('Duration: {}'.format(end_time - start_time))

