#Import necessary libraries
from flask import Flask, render_template, request, send_file

import numpy as np


import os
import scipy
import pandas as pd
import joblib
import chemparse as cp 
import sys
import csv
import datetime
from datetime import date
import time

# Create flask instance
app = Flask(__name__)

#Set Max size of file as 10MB.
# app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024

# #Allow files with extension png, jpg and jpeg
# ALLOWED_EXTENSIONS = ['gz']
# def allowed_file(filename):
#     return '.' in filename and \
#            filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# def init():
#     global graph
#     graph = tf.compat.v1.get_default_graph()





@app.route("/", methods=['GET', 'POST'])
def home():

    return render_template('home.html', input="CH4\nC2H6\nC3H8", result="1.120, 1.744, 2.205")

@app.route("/csvDownload", methods = ['GET', 'POST'])
def csvDownload():

    input = request.form["text"]
    result = request.form["result"]
    result = [x.strip() for x in result.split(',')]
    input = [x.strip() for x in input.split('\n')]
    name = 'csv-storage/predictions' + str(time.time()) + '.csv'
    print(result)
    with open(name, 'w') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(["Input", "Result"])
            for i in range(len(input)):
                print(result[i])
                filewriter.writerow([str(input[i]), str(result[i])])
            path = "/Examples.pdf"
    return send_file(name, as_attachment=True)


@app.route("/about", methods = ['GET','POST'])
def about():
    print(print(request.form["text"]))
    return ""

@app.route("/predict", methods = ['GET','POST'])
def predict():

  if request.method == 'POST':
  #  molecular_formula = request.text['input']
    try:
        input = [x.strip() for x in request.form["text"].split('\n')]
        results = []

        single_compound = 1
        compound_list = 0

        molecular_formula = request.form["text"]
        formula = request.form["text"]

        MFLOGP = joblib.load('MFLOGP.sav')
        scale_X = joblib.load('scale_X.sav')
        scale_y = joblib.load('scale_y.sav')

        
        #'C20H40' #"""Single Compound Only"""#

        #Only required if you are analyzing a list of compounds in an excel document#
        file_dir = r'Input File Path Here'  
        sheetname = 'Input Sheet Name Here'

        #----------------------------------------------------------------------------#
        if (single_compound) | ((single_compound !=1) & (compound_list !=1)) | (single_compound & compound_list):
            for i in range(len(input)):
                if ((single_compound !=1) & (compound_list !=1)):
                    print('No run option chosen, default single compound was run')
                    single_compound = 1
                    
                if (single_compound & compound_list):
                    print('Both run option chosen, default single compound was run')
                    single_compound = 1
                    compound_list = 0

                molecular_formula = pd.DataFrame([input[i]],columns = ['Formula'])
                
                features = molecular_formula['Formula'].apply(cp.parse_formula)
                features = pd.json_normalize(features)
                features = features.fillna(0)
                
                elements = ['C','H','N','O','S','P','F','Cl','Br','I']
                compound = pd.DataFrame(0,index = np.arange(1),columns = [elements])
                
                for ii in range(0,len(features.columns)):
                    if features.columns[ii] == 'C':
                        compound['C'] = features['C']
                
                    elif features.columns[ii] == 'H':
                        compound['H'] = features['H']
                    
                    elif features.columns[ii] == 'N':
                        compound['N'] = features['N']
                        
                    elif features.columns[ii] == 'O':
                        compound['O'] = features['O']
                        
                    elif features.columns[ii] == 'S':
                        compound['S'] = features['S']
                        
                    elif features.columns[ii] == 'P':
                        compound['P'] = features['P']
                        
                    elif features.columns[ii] == 'F':
                        compound['F'] = features['F']
                        
                    elif features.columns[ii] == 'Cl':
                        compound['Cl'] = features['Cl']
                        
                    elif features.columns[ii] == 'Br':
                        compound['Br'] = features['Br']
                        
                    elif features.columns[ii] == 'I':
                        compound['I'] = features['I']
                        
                    else:
                        sys.exit('Incompatible Formula')   
                
                compound_prediction = scale_y.inverse_transform((MFLOGP.predict(pd.DataFrame(scale_X.transform(compound),columns = elements))).reshape(-1,1))
                results.append("{:.3f}".format(compound_prediction[0][0]))

        # elif compound_list:
            
        #     data = pd.read_excel(file_dir,sheet_name = sheetname)
            
        #     features = data['Formula'].apply(cp.parse_formula)
        #     features = pd.json_normalize(features)
        #     features = features.fillna(0)
            
        #     if len(features.columns) > 10:
        #         sys.exit('Incompatible Formula')   
        #     elif len(features.columns)<10:
        #         elements = ['C','H','N','O','S','P','F','Cl','Br','I']
        #         compound = pd.DataFrame(features, columns = elements)
        #         compound = compound.fillna(0)
                
        #     else:
        #         compound = features

        #     MFLOGP = joblib.load('MFLOGP.sav')
        #     scale_X = joblib.load('scale_X.sav')
        #     scale_y = joblib.load('/Users/taruneswar/Documents/Front-end Web Development/NSF CEDAR Lab/web-app/scale_y.sav')

        #     compound_prediction = None
        # try:
        #     MFLOGP = joblib.load('MFLOGP.sav')
        #     scale_X = joblib.load('scale_X.sav')
        #     scale_y = joblib.load('/Users/taruneswar/Documents/Front-end Web Development/NSF CEDAR Lab/web-app/scale_y.sav')
        # except:
        #     print('Please ensure MFLOGP.sav, scale_X.sav, and scale_Y.sav are all in the same directory as MFLOGP_Run_Code.py')
        
        

        # try:

        #     compound_prediction = scale_y.inverse_transform((MFLOGP.predict(pd.DataFrame(scale_X.transform(compound),columns = elements))).reshape(-1,1))
        # except:
        #     print("none")

        # if single_compound:
        #     print(compound_prediction)
        # elif compound_list:
        #     predictions = compound_prediction
        #     print(predictions)
        finalIn = ""
        for i in range(len(input)):
            if(i == len(input)-1):
                finalIn += (str(input[i]))
                break
            else:
                finalIn += (str(input[i]) + "\n")

        finalOut = ""
        for i in range(len(results)):
            if(i == len(results)-1):
                finalOut += (str(results[i]))
                break
            else:
                finalOut += (str(results[i]) + ", ")

        return render_template('home.html', input = finalIn, result= finalOut)
    except:
        return render_template('home.html', input = "", result= "Error: Please enter a valid input")




@app.route("/clear", methods = ['GET','POST'])
def clear():

  return render_template('home.html', input = "", result= "")
        

    




# @app.route("/predict1", methods = ['GET','POST'])
# def predict1():

#     if request.method == 'POST':
#         file = request.files['file']
#         if file and allowed_file(file.filename):
#             filename = file.filename
#             file_path = os.path.join('static/images', filename)
#             file.save(file_path)
#             img = read_image1(file_path)
#             # Predict the class of an image
#            # graph = tf.compat.v1.get_default_graph()
#             with graph.as_default():
#                 model = load_model('NOVEL_schizophrenia_resting_state_dataset_t1_2d.h5')
               
          
#                 img = np.expand_dims(img, axis=0)


#                 predictions = model.predict(img)
#                 print(predictions)
#                 label_index = predictions
#              #   class_prediction = model1.predict_classes(img)
#               #  print(class_prediction)

#             #Map apparel category with the numerical class
#             if label_index == 0.9999999:
#               result = "Healthy TRS MRI Slice"
#             elif label_index != 0.9999999:
#               result = "SZD"
  
#             return render_template('szd.html', result = result)
        

#     return render_template('szd.html')



if __name__ == "__main__":
   
    app.run(debug=True)

