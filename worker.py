#find newly added job and do the job
#remove old item and remove
import numpy as np
import os, time, sys, code, time, shutil, pickle,textdistance
from stat import S_ISREG, ST_CTIME, ST_MODE
import RNA
import time
import docker
import zipfile

import mimetypes, smtplib
from email.mime.text import MIMEText

dir_path = './static/results/'
data_path = "./dataset/"
# get all entries in the directory

y_mapping = {"CAS-TypeIE":0,"CAS-TypeIIC":1,"CAS-TypeIB":2,"CAS-TypeIC":3,"CAS-TypeIF":4,"CAS-TypeIIIB":5,"CAS-TypeIA":6,"CAS-TypeIIA":7,"CAS-TypeIIIA":8,"CAS-TypeID":9,"CAS-TypeIIID":10,"CAS":11,"CAS-TypeIU":12,"CAS-TypeVIB1":13,"CAS-TypeIIB":14,"CAS-TypeVA":15,"CAS-TypeVIA":16,"CAS-TypeIIIC":17,"CAS-TypeIV":18,"CAS-TypeVB":19,"CAS-TypeVIC":20,"CAS-TypeVIB2":21}
y_label = ["CAS-TypeIE","CAS-TypeIIC","CAS-TypeIB","CAS-TypeIC","CAS-TypeIF","CAS-TypeIIIB","CAS-TypeIA","CAS-TypeIIA","CAS-TypeIIIA","CAS-TypeID","CAS-TypeIIID","CAS","CAS-TypeIU","CAS-TypeVIB1","CAS-TypeIIB","CAS-TypeVA","CAS-TypeVIA","CAS-TypeIIIC","CAS-TypeIV","CAS-TypeVB","CAS-TypeVIC","CAS-TypeVIB2"]

print("Worker is getting started...")

while(True):
    results = os.listdir(dir_path)
    for result in results:
        path = dir_path + result
        #time_created = os.stat(path).st_ctime
        #if (time.time() - time_created)/3600/24/15 >1: #remove results
        #    shutil.rmtree(path)
        current_pth = os.getcwd()
        try:
            if os.path.exists(path):
                if len(os.listdir(path)) == 4:
                    st_time = time.time()
                    # Check if protein file is available
                    if 'protein.fasta' in os.listdir(path) and 'alphafold.pdb' not in os.listdir(path):
                        command = f"python /home/gugcweb/alphafold/docker/run_docker.py \
                                --fasta_paths=/home/gugcweb/Cas-CRISPR\ matcher\ server/static/results/{result}/protein.fasta \
                                --max_template_date=1950-05-14 \
                                --model_preset=monomer \
                                --db_preset=reduced_dbs \
                                --output_dir=/home/gugcweb/Cas-CRISPR\ matcher\ server/static/results/{result}/ \
                                --data_dir=/media/new_partition/data_alphafold/"
                        os.system(command)
                        os.system(f'mv {path}/protein/ranked_0.pdb {path}/alphafold.pdb')
                    else:
                        pass
                    print("alphafold finished", time.time()- st_time)
                    # Check if RNA sequence or RNA pdb file
                    if 'ternary.pdb' in os.listdir(path):
                        current_pth = os.getcwd()
                        os.chdir(path)

                        # Remove (incomplete) first residue atoms 
                        with open("ternary.pdb", "r") as f:
                            lines = f.readlines()
                        res_num_first = lines[1].split()[5]
                        lines_refined = [i for i in lines[1:-2] if i.split()[5] != res_num_first]

                        with open("ternary_refined.pdb", "w") as f:
                            f.writelines(lines_refined)
                    else:
                        print(path)
                        primary_seq = open(path + "/" + "primary.fasta").readlines()[-1].strip()

                        secondary_seq = RNA.fold(primary_seq)
                        f_out = open(path + "/secondary.secstruct","w")
                        f_out.write(secondary_seq[0].lower() + "\n" + primary_seq.lower() + "\n")
                        f_out.close()

                        print("Start to build ternary structure")

                        current_pth = os.getcwd()
                        os.chdir(path)

                        #Build ternary structure
                        command = "/media/ext_usb_drive/dataset_lib/rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/rna_denovo.static.linuxgccrelease -fasta primary.fasta -secstruct_file secondary.secstruct -out:file:silent intermediate.out -minimize_rna true"
                        print(command)
                        os.system(command)

                        command2 = """grep "ˆSCORE:" intermediate.out | grep -v description | sort -nk2 | head -n 500 | sort -nk24 | head -n 1 | awk '{print $NF}'"""
                        print(command2)
                        os.system(command2)
                        
                        command3 = f"""/media/ext_usb_drive/dataset_lib/rosetta/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/extract_pdbs.static.linuxgccrelease -in::file::silent intermediate.out -tags $TAG"""
                        print(command3)
                        os.system(command3)

                        #remove first element(neucleotide) of ternary.pdb 
                        
                        f_out = open("ternary.pdc","w")
                        for line in open("S_000001.pdb").readlines():
                            if line[:4] == "ATOM":
                                if int(line[22:26]) == 1:
                                    continue
                            f_out.write(line)

                        os.system("rm *.out")
                        os.system("mv S_000001.pdb ternary.pdb")
                        os.system("mv ternary.pdc ternary_refined.pdb")

                    # Get Protein PDB file names
                    pdbs = [i for i in os.listdir('.') if i.endswith('.pdb')]
                    print(pdbs)
                    protein_pdb = [i for i in pdbs if i != 'ternary.pdb' and i !='ternary_refined.pdb'][0]
                    
                    
                    command4 = f"""../../../HDOCKlite/hdock {protein_pdb} ternary_refined.pdb -out Hdock.out"""
                    print(command4)
                    os.system(command4)
                    command5 = f"""../../../HDOCKlite/createpl Hdock.out top1.pdb -nmax 10 -complex -models"""
                    os.system(command5)
                    
                    docking_zip = zipfile.ZipFile('docking_results.zip', 'w')
                    model_list = ['model_1.pdb', 'model_2.pdb', 'model_3.pdb', 'model_4.pdb', 'model_5.pdb', 'model_6.pdb', 'model_7.pdb', 'model_8.pdb', 'model_9.pdb', 'model_10.pdb']
                    for i in model_list:
                        docking_zip.write(i, compress_type=zipfile.ZIP_DEFLATED)
                    docking_zip.close()

                    task_file = open("task_name.txt","r")
                    task_name = task_file.readlines()[0]
                    task_file.close()

                    email_file = open('email.txt', 'r')
                    user_email = email_file.readlines()[0]
                    email_file.close()

                    smtp = smtplib.SMTP('smtp.gmail.com', 587)
                    smtp.starttls()
                    smtp.login('crisprcasdocker@gmail.com', 'password')
                    msg_html = f"""\
                    <html>
                    <head></head>
                        <body>
                            <p><b>Your task is finished!</b></p>
                            <p>Email: {user_email}</p>
                            <p>Task Name: {task_name}</p>
                            <p>Use the link below to get your result.</p>
                            <p>http://www.crisprcasdocker.org/results/{result}</p>
                        </body>
                    </html>
                    """
                    msg = MIMEText(msg_html, 'html')
                    msg['To'] = user_email
                    msg['Subject'] = 'CCD - Your Job is Finished!' 
                    smtp.sendmail('crisprcasdocker@gmail.com', user_email, msg.as_string())
                    smtp.quit()

                    os.chdir(current_pth)
        except:
            error_file = open('error.log', 'w')
            error_file.write('Something went wrong')
            error_file.close()
            os.chdir(current_pth)
            
    time.sleep(5)
    print('sleep')
    #code.interact(local=dict(globals(), **locals()))

    #../source/bin/rnadenovo.static.linuxgccrelease -fasta chunk002_1lnt.fasta -secstructfile chunk002_1lnt.secstruct -nstruct 2 -out:file:silent test.out -minimize_rna -dump=True
    #../source/bin/rnadenovo.static.linuxgccrelease

