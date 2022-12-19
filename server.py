from flask import Flask,render_template,request, redirect, send_file
import random, code, hashlib, os,smtplib, dash, flask,json
from email.mime.text import MIMEText
from multiprocessing import Process
import plotly
import plotly.express as px
import pandas as pd
import numpy as np
import textdistance, pickle
import plotly.graph_objects as go
#import dash_html_components as html
#import dash_core_components as dcc
# from dash import dcc
import urllib.request as urlreq
#from dash import Dash
#import dash_bio as dashbio
#import dash_bootstrap_components as dbc
#from dash_bio.utils import PdbParser, create_mol3d_style
from PIL import Image
from random_word import RandomWords
from Bio import SeqIO

app = Flask(__name__)
external_scripts=[{'src': "https://3Dmol.org/build/3Dmol-min.js"}]
#app_dash = dash.Dash(__name__, server=app, url_base_pathname='/show_dashboard/')
#app_dash.config.suppress_callback_exceptions = True
#app_dash.layout = html.Div('Wrong direction')

result_path = "/home/gugcweb/Cas-CRISPR matcher server/static/results/"

def wait_template(query_id):
	return f'''
	<!doctype html>
	<html>
		<head>
			<META HTTP-EQUIV="refresh" CONTENT="30">
		</head>
		<body>
			<h1>Your task {query_id} is now working</a></h1>

			Please wait while the job is running. This page keeps refreshing automatically every 30 seconds.</br>
			A link to the result page will be sent to your email when it is done.

			<ul>
			<li>Go back to <a href="/">Main page</a></li>
			</ul>
		</body>
	</html>

	'''

@app.route("/results", methods = ["GET","POST"])
def result():
	if request.method == "GET":
		return render_template('result.html')
	else:
		email = request.form['mail']
		task_name = request.form['task_name']
		email_hash = hashlib.sha256(email.encode()).hexdigest()[:10]
		task_hash = hashlib.sha256(task_name.encode()).hexdigest()[:10]
		
		query_id = email_hash + task_hash
		url = '/results/'+ str(query_id)

		return redirect(url)


@app.route("/contact")
def contact():
	return render_template('contact.html')


@app.route("/help")
def help():
	return render_template('help.html')


@app.route("/", methods = ["GET","POST"])
def index():
	if request.method == "GET": # with sequence checking
		return render_template('main.html')
	else:
		global result_path

		email = request.form['mail']
		
		if request.form['task_name']:
			task_name = request.form['task_name']
		else:
			r = RandomWords()
			task_name = r.get_random_word() # Return a single random word

		email_hash = hashlib.sha256(email.encode()).hexdigest()[:10]
		task_hash = hashlib.sha256(task_name.encode()).hexdigest()[:10]
		
		query_id = email_hash + task_hash
		path = result_path+query_id

		#save seq 
		if not os.path.exists(path):
			os.mkdir(path)

		f_out = open(os.path.join(path,"task_name.txt"),"w")
		f_out.write(task_name)
		f_out.close()

		f_out = open(os.path.join(path,"email.txt"),"w")
		f_out.write(email)
		f_out.close()


		# Get RNA sequence or PDB file
		rna_sequence = request.form.get('rna_seq')
		rna_file_fasta = request.files.get('rna_file_fasta')
		rna_file_pdb = request.files.get('rna_file_pdb')

		print(rna_sequence, rna_file_fasta, rna_file_pdb)
		if rna_file_pdb:
			rna_file_pdb.save(os.path.join(path, 'ternary.pdb'))
		elif rna_file_fasta:
			fasta_sequence = rna_file_fasta.read().decode('utf-8')
			print(fasta_sequence)
			f_out = open(os.path.join(path,"primary.fasta"),"w")
			f_out.write(fasta_sequence.lower().replace("t","u").replace('*', ''))
			f_out.close()

		elif not rna_file_pdb and not rna_file_fasta and rna_sequence:
			f_out = open(os.path.join(path,"primary.fasta"),"w")
			f_out.write(rna_sequence.lower().replace("t","u").replace('*', ''))
			f_out.close()
		else:
			print("No RNA data is uploaded via form!!!")  # The query will not be processed in this case

		# get PDB file or ID
		cas_id = request.form['cas_id']
		cas_file_pdb = request.files['cas_file_pdb']
		cas_file_fasta = request.files['cas_file_fasta']
		is_caspdb = request.form.get('caspdb')

		if is_caspdb:
		
			if cas_file_pdb and cas_id:
				cas_name = cas_file_pdb.filename
				cas_file_pdb.save(os.path.join(path, cas_name))
			elif cas_file_pdb and not cas_id:
				cas_name = cas_file_pdb.filename
				cas_file_pdb.save(os.path.join(path, cas_name))
			elif not cas_file_pdb and cas_id:
				path_temp = path.replace(' ', '\ ')
				command = f"""wget -P {path_temp} -N https://files.rcsb.org/download/{cas_id}.pdb"""
				os.system(command)
			else:
				print("No protein data is uploaded via form!!!") # The query will not be processed in this case

		else:
			if cas_file_fasta:
				cas_name = cas_file_fasta.filename
				cas_file_fasta.save(os.path.join(path, cas_name))
			else:
				protein_seq = request.form['prot_seq']
				f_out = open(os.path.join(path,"protein.fasta"),"w")
				f_out.write('>query\n')
				f_out.write(protein_seq.replace('*', ''))
				f_out.close()

		url = '/results/'+ str(query_id)
		return redirect(url)

@app.route("/results/<query_id>",methods = ["GET"])
def get_results(query_id):
	#check whether the quesy is done
	global result_path

	if query_id in os.listdir(result_path):
		path = result_path + query_id
		item_list = os.listdir(path)
		print(item_list)
		if item_list != None and len(item_list) > 3 and 'model_1.pdb' in item_list:
			#show html page for everything
			print("there is somthing in the result")
			cas = [i for i in item_list if i.endswith('.pdb')]
			print(cas)
			cas.remove('ternary.pdb')
			cas.remove('ternary_refined.pdb')
			for i in range(1,11):
				cas.remove(f'model_{i}.pdb')
			cas_base_path =  path + '/'
			cas_path = os.path.join(path, cas[0])
			rna_path = os.path.join(path, 'ternary_refined.pdb')
			model_path = os.path.join(path, 'model_1.pdb')

			dir_path = '/home/gugcweb/Cas-CRISPR matcher server/static/results/'
			data_path = "/home/gugcweb/Cas-CRISPR matcher server/dataset/"
			rna_seq_path = os.path.join(dir_path, query_id)

			# Graph 1 
			y_mapping = {"CAS-TypeIE":0,"CAS-TypeIIC":1,"CAS-TypeIB":2,"CAS-TypeIC":3,"CAS-TypeIF":4,"CAS-TypeIIIB":5,"CAS-TypeIA":6,"CAS-TypeIIA":7,"CAS-TypeIIIA":8,"CAS-TypeID":9,"CAS-TypeIIID":10,"CAS":11,"CAS-TypeIU":12,"CAS-TypeVIB1":13,"CAS-TypeIIB":14,"CAS-TypeVA":15,"CAS-TypeVIA":16,"CAS-TypeIIIC":17,"CAS-TypeIV":18,"CAS-TypeVB":19,"CAS-TypeVIC":20,"CAS-TypeVIB2":21}
			y_label = ["CAS-TypeIE","CAS-TypeIIC","CAS-TypeIB","CAS-TypeIC","CAS-TypeIF","CAS-TypeIIIB","CAS-TypeIA","CAS-TypeIIA","CAS-TypeIIIA","CAS-TypeID","CAS-TypeIIID","CAS","CAS-TypeIU","CAS-TypeVIB1","CAS-TypeIIB","CAS-TypeVA","CAS-TypeVIA","CAS-TypeIIIC","CAS-TypeIV","CAS-TypeVB","CAS-TypeVIC","CAS-TypeVIB2"]
			y_label_no_cas = ["IE","IIC","IB","IC","IF","IIIB","IA","IIA","IIIA","ID","IIID","CAS","IU","VIB1","IIB","VA","VIA","IIIC","IV","VB","VIC","VIB2"]
			
			print(1.1)
			print(os.getcwd())
			with open(os.path.join(rna_seq_path, 'primary.fasta')) as f:
				rna_sequence = f.readlines()[0].upper().rstrip('\n')
			target_seq = rna_sequence
			print(target_seq)
			

			print(1.11)
			lines = open(data_path +"eigen_seq.csv").readlines()
			dists = [] 
			labels = []
			seqs = []
			print(1.12)
			for idx, line in enumerate(lines):
				#print(idx, len(lines))
				parse = line.split(",")
				name, cas_type, occurances, seq = line.split(",") 
				labels.append(cas_type)
				dists.append(textdistance.hamming(target_seq, seq.strip()))
				seqs.append(seq)

			# rearrange [([line], dists)] according to distance in ascending order
			res = sorted(zip(lines, dists), key=lambda x:x[1])  # res?	
			print(1.2)

			#Generating fasta format for MSA
			ret_fasta = ">sp|id|target\n" + target_seq +"\n"  # ret?

			for item, dist in res[:10]:
				name, cas_type, occurances, seq = item.split(",")
				ret_fasta +=">sp|id|" + name +"_" + cas_type.replace("CAS-Type","") + "_" + str(dist) + "\n" + seq.strip() + "\n"
			
			#Load tsne for scatter plot 
			#-> need to change or update? like to tsne everytime?
			#-> This code cannot detect the worst sequence case. like too long seq too weired case which min dist is larger than 10
			X2 = pickle.load(open(data_path+"hamming_eigen_tsne.pickle","rb"))
			hue = [y_mapping[item] for item in labels]
			hue = [item.replace("CAS-Type","") for item in labels]

			dist_threshold = sum(sorted(dists)[:10])/10

			pred_location = X2[np.where(np.asarray(dists) < dist_threshold)].mean(axis  = 0)
			pred_label = 'Target_' + labels[np.argmin(dists)].replace("CAS-Type","")
			print(1.3)
			if np.isnan(pred_location[0]): #In the case when there is no less than distance 5
				pred_location = np.asarray([-100,-100])
				pred_label = "Too far from existing sequences"
			
			df = pd.DataFrame({"t-SNE 1":X2[:,0],"t-SNE 2":X2[:,1],"label":hue,"seq":seqs })

			fig = go.Figure()
			#df.to_csv("temp.csv")
			fig = px.scatter(df, x="t-SNE 1", y="t-SNE 2", color= "label", symbol = "label",hover_data=["seq"])
			fig.add_trace(go.Scatter(x = [pred_location[0]], y = [pred_location[1]],		# Edited
								mode = 'markers',
								marker = dict(color='red', size = 15, symbol ="star-dot"),
								name = pred_label))
			fig.update_layout(xaxis_range=[-50+pred_location[0],pred_location[0]+50],yaxis_range=[-50+pred_location[1],50+pred_location[1]])
			graph1 = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


			# Graph 2
			neighbor_dist = np.zeros((22,)) + 0.1
			#neighbor_dist = {"IE":0,"IIC":0,"IB":0,"IC":0,"IF":0,"IIIB":0,"IA":0,"IIA":0,"IIIA":0,"ID":0,"IIID":0,"CAS":0,"IU":0,"VIB1":0,"IIB":0,"VA":0,"VIA":0,"IIIC":0,"IV":0,"VB":0,"VIC":0,"VIB2":0}
			for line, dist in res[:30]:
				cl = line.split(",")[1]
				neighbor_dist[y_mapping[cl]] += 1
			#fig = px.bar(data_canada, x='year', y='pop')
			neighbor_dist = neighbor_dist/sum(neighbor_dist)
			fig_bar = px.bar(y=neighbor_dist, x=y_label_no_cas,color= y_label_no_cas)
			
			fig_bar.update_layout(xaxis_title="Class name", yaxis_title="Maximum Ratio of Class")
			graph2 = json.dumps(fig_bar, cls=plotly.utils.PlotlyJSONEncoder)

			path = '/static/results/' + query_id
			cas_base_path =  path + '/'
			cas_path = os.path.join(path, cas[0])
			rna_path = os.path.join(path, 'ternary_refined.pdb')
			model_path = os.path.join(path, 'model_1.pdb')


			return render_template('plot_result.html', query_id = query_id, cas_path=cas_path, cas_base_path=cas_base_path, rna_path=rna_path, model_path=model_path, graphJSON1 = graph1, graphJSON2 = graph2)


		else:
			#show wating page
			f_out = open(os.path.join(path,"task_name.txt"),"r")
			task_name = f_out.readlines()[0]
			f_out.close()
			return wait_template(task_name)
	else:
		#wrong direction.
		return "Wrong direction" # TO DO: make independent HTML template


@app.route("/download",methods = ["POST"])
def model_download():
	file_path = request.form.get('file_path').strip('/')
	print(file_path)
	if '.zip' in file_path:
		return send_file(file_path,
					mimetype='application/zip',
					attachment_filename='docking_results.zip',
					as_attachment=True)
	else:
		if 'model' in file_path:
			return send_file(file_path,
						mimetype='chemical/pdb',
						attachment_filename='docking_result.pdb',
						as_attachment=True)
		elif 'ternary' in file_path:
			return send_file(file_path,
						mimetype='chemical/pdb',
						attachment_filename='crRNA_3D.pdb',
						as_attachment=True)
		else:
			return send_file(file_path,
						mimetype='chemical/pdb',
						attachment_filename='protein_3D.pdb',
						as_attachment=True)



if __name__ == '__main__':
	app.run(host = '0.0.0.0', port=8888)#debug -> 변경된 부분이 생기면 자동으로 껏다켜줌 나중에 삭제 필요 
