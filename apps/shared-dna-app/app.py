import os, re, glob
import pandas as pd
import numpy as np
from flask import Flask, flash, request, render_template, redirect, url_for, session
from werkzeug.utils import secure_filename
from functions import *
global configs
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

app = Flask(__name__, static_url_path='/output')

UPLOAD_FOLDER = 'genomes'
RESULT_FOLDER = 'output'

app._static_folder = RESULT_FOLDER

app.secret_key = "secret key"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['RESULT_FOLDER'] = RESULT_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 2000 * 1024 * 1024
app.config['STATIC_FOLDER'] = RESULT_FOLDER

ALLOWED_EXTENSIONS = set(['txt', 'zip', 'gz'])
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def upload_form():
    return render_template('./upload.html')

@app.route('/', methods=['POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        file1 = None
        if 'genome1' in request.files:
            file1 = request.files['genome1']
        if file1:
            if (not allowed_file(file1.filename)) :
                flash('The genome file is not valid. The allowed file types are txt, zip, or gz')
                return redirect(request.url)
            filename1 = secure_filename(file1.filename)
            file1.save(os.path.join(app.config['UPLOAD_FOLDER'], filename1))
            flash('Files successfully uploaded')
            genome1 = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
            name1 = request.form['name1']
            if name1 == '':
                name1 = 'Me'
            flash("Now checking all the traits for " + name1 + ' ...')
            session['name1'] = name1
            session['genome1'] = genome1
            configs['data_type'] = request.form['format'] 
            configs['background'] = './resources/' + request.form['background'] + '_distributions.tsv.gz' 
            configs['models'] = './resources/gwas_filtered' + request.form['models'] + '.csv.gz' 
            # remove old files
            for oldfile in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], name1 + '*')):
                os.remove(oldfile)
            return redirect('/analyzeMe')

        else:
            flash('At least one valid genome file is required. The allowed file types are txt, zip, or gz')
            return redirect(request.url)

@app.route('/analyzeMe')
def process_data():
    if configs['data_type'] == 'scores':
        session['matched'] = session['genome1']
    else:
        name1 = session['name1']
        genome1 = session['genome1']
        matched  = matchOne(genome1, name1)
        matched_file = os.path.join(app.config['RESULT_FOLDER'], name1 + "_matched.tsv.gz")
        matched.to_csv(matched_file, sep='\t', compression='gzip')
        session['matched'] = matched_file 
        os.remove(genome1)
    flash('Analysis completed! Now you can use keywords to explore yourself.')
    return render_template('show_result.html', tables=[], png='')



@app.route('/show_result', methods=['POST'])
def search_result():
    if request.method == 'POST':
        term = request.form['term']
    if term == '':
        term = 'Brain'
    matched = session['matched']
    matched = pd.read_csv(matched, sep='\t', index_col=0)
    for t in term.split(' '):
        matched = matched.filter(regex=re.compile(t, re.IGNORECASE), axis=0)
    name1 = session['name1']
    if matched.shape[0] > 10:
        matched = matched.iloc[0:10,:]
        flash('Too many results for ' + term + ', only showing 10. Try to use multiple keywords for searching.')
    if matched.shape[0] == 0:
        flash('No results for ' + term + '. Try sth. different.')
        return render_template('show_result.html', tables=[], png='')
    fig = create_figure(matched)
    png_file = os.path.join(app.config['RESULT_FOLDER'], name1 + '_' + term.replace(' ', '_') + ".png")
    fig.savefig(png_file, bbox_inches ='tight' )
    return render_template('show_result.html', tables=format_match_table(matched.iloc[:,0:5]), png=png_file)

def create_figure(matched):
    """
    create a boxplot and show user data
    """
    data = matched[['mean', '50%', '5%', '25%', '75%', '95%']]
    data = data.rename(columns={'50%':'median', '5%':5, '25%':25, '75%':75, '95%':95})
    data['user'] = matched.iloc[:,0]
    data = data.div(matched['max'], axis=0)
    return boxplot(data)

@app.route('/download')
def download_data():
    matched = session['name1'] + "_matched.tsv.gz"
    return render_template('download.html', results=matched)

if __name__ == '__main__':
    app.run(host="0.0.0.0")
