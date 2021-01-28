#!/usr/bin/env python3

"""
A very simple HTTP server in python for logging requests
Usage::
    ./server.py [<port>]
"""

from http.server import BaseHTTPRequestHandler, HTTPServer
import logging
import cgi
import urllib.request
import tarfile
import tempfile
import smtplib
import subprocess
import glob
import operator
from concurrent.futures import ThreadPoolExecutor
from numpy import *

from email.message import EmailMessage
from email.headerregistry import Address


GROUP_NAMES = ['SBROD-server','SBROD-plus']
MY_EMAIL = "xxx@gmail.com"
GROUP_NAME = MY_EMAIL

QA_METHODS = ['./model1/sbrod', './model2/sbrod']
AUTHORS = ['xxxx-xxxx-xxxx', 'yyyy-yyyy-yyyy']
METHODS = [
    'SBROD (trained on servers from CASP5-5-12) + NMA models',
    'SBROD (trained on all predictions from CASP5-5-12) + NMA models'
]

SENDER_EMAIL = "xxxxx.xx@gmail.com"
SENDER_PASS = "XXXXXXXXXXXX"


def send_email(From, To, Subject, Content):
    # Create the base text message
    msg = EmailMessage()
    msg['Subject'] = Subject
    msg['From'] = From
    msg['To'] = To
    msg['CC'] = MY_EMAIL
    msg.set_content(Content)

    logging.info('Email to send:%s\n', msg)

    # Send the message via gmail
    with smtplib.SMTP('smtp.gmail.com', 587) as server:
        server.ehlo()
        server.starttls()
        server.login(SENDER_EMAIL, SENDER_PASS)
        server.send_message(msg)


def prepare_output_from_scores(text, method):
    target = text.split()[0].split('/')[-2]
    ar=[]
    scores = {}
    for a, b in zip(text.split()[0::2], text.split()[1::2]):
        scores[a.split('/')[-1]] = float(b)
        ar.append(float(b))

    ar = array(ar)

    model = 1
    if(next(iter(scores))[:6] == 'server'):
        model = 1
    else:
        model = 2

    output = '''PFRMAT QA
TARGET '''+target+'''
AUTHOR '''+AUTHORS[method]+'''
METHOD '''+METHODS[method]+'''
MODEL '''+str(model)+'''
QMODE 1
'''
    for key, value in sorted(scores.items(), key=operator.itemgetter(1), reverse=1):
        output += '%s %f\n'%(key, value)

    output += 'END\n'
    return output


def process(target, tarball, email, model):
    LOG_MESSAGE = "Request ({}, {}, {}, {}): %s\n".format(target, tarball, email, model)
    logging.info(LOG_MESSAGE, "start processing request")

    try:
        # Download the file
        ftpstream = urllib.request.urlopen(tarball, timeout=10)
        loadedtarfile = tarfile.open(fileobj=ftpstream, mode="r|gz")

        with tempfile.TemporaryDirectory() as tmpdirname:
            loadedtarfile.extractall(tmpdirname)

            # Score all proteins in `tmpdirname`
            pdb_files = glob.glob(tmpdirname + '/*/*')
            logging.info(LOG_MESSAGE, "start scoring files: {}".format(pdb_files))

            if (model=='servers'):
            	QA_METHOD = QA_METHODS[0]
            	method = 0
            else :
            	QA_METHOD = QA_METHODS[1]
            	method = 1

            try:
            	result = subprocess.check_output([QA_METHOD] + pdb_files + ['--scale'])
            except subprocess.CalledProcessError as resultEx:   
    	        raise RuntimeError("ERROR: protein QA exit code: {} {}".format(
                                        resultEx.returncode, resultEx.grepexc.output))

        scores = result.decode('utf-8')

        output = prepare_output_from_scores(scores, method)
        logging.info(LOG_MESSAGE, 'Content to send:\n' + output)

        send_email(MY_EMAIL, email, target + ' scores by ' + GROUP_NAME, output)
    except Exception as e:
        logging.exception(LOG_MESSAGE, e)


class S(BaseHTTPRequestHandler):

    def do_POST(self):
        form = cgi.FieldStorage(fp=self.rfile,
                                headers=self.headers,
                                environ={'REQUEST_METHOD': 'POST'})

        logging.info("POST request,\nPath: %s\nHeaders:\n%s\n\nBody:\n%s\n",
                     str(self.path), str(self.headers), form)

        try:
            target = form['target'].value
            tarball = form['tarball'].value
            email = form['email'].value
            model = form['model'].value
            self.send_response(200)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write("{} - query received by {}\n".format(target, GROUP_NAME).encode('utf-8'))
            self.executor.submit(process, target, tarball, email, model)
        except Exception as e:
            logging.exception(str(e))
            self.send_response(400)
            self.send_header('Content-type', 'text/html')
            self.end_headers()
            self.wfile.write("A bad query received by {}\n".format(GROUP_NAME).encode('utf-8'))


def run(port=8080):
    executor = ThreadPoolExecutor(1)
    S.executor = executor

    logging.basicConfig(level=logging.INFO)
    server_address = ('', port)
    httpd = HTTPServer(server_address, S)
    logging.info('Starting httpd on port {}...\n'.format(port))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    httpd.server_close()
    logging.info('Stopping httpd...\n')


if __name__ == '__main__':
    from sys import argv

    if len(argv) == 2:
        run(port=int(argv[1]))
    else:
        run()
