
from flask import Flask
from flask import render_template
from werkzeug.wrappers import Request, Response
import pandas as pd
import json


app = Flask(__name__)

@app.route('/')
def gcg():

    url = 'plotdata.pkl'
    plotdata = pd.read_pickle(url)

    title = "Time Distribution"
    #temps['Month'] = temps['Month'].astype(str)
    #temps['Avg'] = (temps['Tmax']+temps['Tmin'])/2

    idx = plotdata.index.tolist()
    d = plotdata.values.tolist()
    c = ["(INSTANCE, DECOMP)"] + plotdata.columns.tolist()
    d.insert(0,c)
    for i in range(len(d)-1):
        inst = idx[i][0].split("/")[-1]
        dec = str(idx[i][1]).split("/")[-1]
        d[i+1] = [f"{inst}; {dec}"] + d[i+1]
    tempdata = json.dumps({'title':title,'data':d})

    return render_template('gcg.html', tempdata=tempdata)

if __name__ == '__main__':
    from werkzeug.serving import run_simple
    run_simple('localhost', 5000, app)
