#!/usr/bin/env python
'''
ArcGet REST
@author: tokeefe
@date: 08.31.2011
@version: 0.2
'''
import time
import io
from datetime import datetime
import re
import os
import sys
import json
import gzip
import shutil
import zipfile
import getpass
import argparse
import traceback
from lxml import etree
import urllib2
from urllib2 import URLError,HTTPError

global progname
progname = os.path.basename(sys.argv[0])
global verbosity
verbosity = 1
global global_t0
global_t0 = time.time()


class ArcGet:
    def __init__(self, args):
        self.args = args
        self.t0 = time.time()
        self.nons = re.compile('.*/')
        self.session_cols = [
            'label',
            'ID',
            'date',
            'time',
            'acquisition_site',
            'scanner',
            'operator',
            'dcmAccessionNumber',
            'dcmPatientId',
            'dcmPatientName',
            'visit',
            'coil',
            'marker',
            'stabilization',
            'consoleversion',
            'xnat:mrSessionData/subject_ID'
            ]
        self.session_cols_lc = []
        for col in self.session_cols:
            self.session_cols_lc.append(col.lower())
        self.passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
        self.auth_handler = urllib2.HTTPBasicAuthHandler()
        self.auth_handler.add_password("XNAT Protected Area",args.host+"/",self.args.username,self.args.password)
        self.opener = urllib2.build_opener(self.auth_handler)
        urllib2.install_opener(self.opener)
    def getArgs(self):
        return self.args
    def getage(self,dob_string,scan_date_string):
        scan_date = datetime.strptime(scan_date_string,"%Y-%m-%d")
        dob_date ="";
        try:
            dob_date = datetime.strptime(dob_string,"%Y-%m-%d")
        except ValueError:
            try:
                dob_date = datetime.strptime(dob_string+"-07-07","%Y-%m-%d")
            except ValueError:
                warn("Unable to determine Age for DOB=\""+dob_date+"\"", verbosity_level = 1)
                return 0
        age = scan_date - dob_date;
        return int(age.days / 365.25)
    def tock(self,t0):
        t1 = time.time()
        info("Time: "+ str(t1 - t0)+" seconds; Total Running time: "+ str(t1 - self.t0)+" seconds", verbosity_level = 2)
    def getRequest(self,url):
        url = self.args.host + url
        info("Getting: "+url, verbosity_level = 3)
        request = None
        response = None
        content = None
        urllib2.install_opener(self.opener)
        request = urllib2.Request(url);
        try:
            request.t0 = time.time()
            response = urllib2.urlopen(request);
        except HTTPError, e:
            self.tock(request.t0)
            if(401 == e.getcode()):
              error("Server said \""+str(e.message)+" Permission Denied\". Do you have access to the data? Are you sure you used the correct username and password?")
            if request is not None:
                error_noexit("HTTP Request Headers: "+str(request.headers))
            else:
                error_noexit("HTTP Request Headers: None")
            if response is not None:
                error_noexit("HTTP Response Headers: "+str(response.headers))
            else:
                error_noexit("HTTP Response Headers: None")
            if content is not None:
                error_noexit("HTTP Response Content: "+str(content))
            else:
                error_noexit("HTTP Response Content: None")
            error("Failed to connect to '"+url+"' "+str(e)+"\nRequest: "+str(request.headers))
        except URLError, e:
            self.tock(request.t0)
            error("Failed to connect to '"+url+"' Reason: "+str(e.reason))
        except Exception, e:
            self.tock(request.t0)
            error("Unknown Error. Can't continue because I'm not smart enough to know what to do.\nError: "+str(e))
        if(401 == response.code):
            self.tock(request.t0)
            error("Server said \""+str(response.code)+" Permission Denied\". Do you have access to the data? Are you sure you used the correct username and password?")
        if(200 != response.code):
            self.tock(request.t0)
            error("Server returned an error (HTTP " + str(response.code) + ") upon attempting to retrieve '"+url+"', don't know what to do about it.")
        response.t0 = time.time()
        return (request,response)
    def geturl(self,url):
        (request,response) = self.getRequest(url)
        content = None
        try:
            content = response.read();
            self.tock(response.t0)
        except Exception, e:
            self.tock(response.t0)
            error("Unknown Error. Can't continue because I'm not smart enough to know what to do.\nError: "+str(e))
        return (response, content)
    def getSessionInfoByLabel(self,label):
        ## --- get Accession ID for MR Session label
        info("Getting Internal ID for MR Session '" + label + "'", verbosity_level = 2)
        url = ("/data/archive/experiments?columns="+",".join(self.session_cols)+"&label=" +
               urllib2.quote(label) + "&format=json")
        # print "self='%s'" % str(self)
        # print "url='%s'" % str(url)
        sessions = self._sessionsFromJsonUrl(url)
        if(len(sessions) > 1):
            error("Multiple sessions returned for MR Session '" + label + "'")
        return sessions
    def getSessionInfoByLabels(self,labels):
        if(len(labels) > 1):
            all_sessions = self.getAllSessionInfo()
            label_dict = {}
            sessions_by_labels = {}
            sessions_by_id = {}
            for label in labels:
                label_dict[label] = True
            for session in all_sessions:
                if session['label'] in label_dict:
                    # If we've seen this label before, and it's not just a duplicate
                    if (session['label'] in sessions_by_labels):
                        if not (session['ID'] in sessions_by_id):
                            # If this is the first time we've run into this, then
                            # print the previously saved session
                            if sessions_by_labels[session['label']] is not False:
                                warn("First session '"+session['label']+"' has ID='"+sessions_by_labels[session['label']]['ID']+"'")
                                # Don't output the first session again.
                                sessions_by_labels[session['label']] = False
                                # Do output this one.
                                warn("Another session '"+session['label']+"' found with ID='"+session['ID']+"'. Adding any way.")
                    else:
                        sessions_by_labels[session['label']] = session
                    sessions_by_id[session['ID']] = session
            sessions = []
            for label in label_dict:
                if not(label in sessions_by_labels):
                    warn("Unable to find MR Session '"+label+"'")
            for session_id in sessions_by_id:
                sessions.append(sessions_by_id[session_id])
            return sessions
        else:
            return self.getSessionInfoByLabel(labels[0])
    def getShowAllSessionInfo(self):
        all_sessions = self.getAllSessionInfo()
        col_sizes = {}
        my_cols = ['Label','ID','Date','Time','Visit' ]
        my_cols_lc = []
        for col in my_cols:
            my_cols_lc.append(col.lower())
        for session in all_sessions:
            for idx in my_cols_lc:
                if idx not in col_sizes:
                    col_sizes[idx] = 5
                if len(session[idx]) > col_sizes[idx]:
                    col_sizes[idx] = len(session[idx])
        print_str = "%-"+ "s %-".join(str(col_sizes[idx]) for idx in my_cols_lc) + "s";
        print print_str % tuple(my_cols);
        for session in all_sessions:
            vals = []
            for col in my_cols_lc:
                # print session
                val = session[col]
                vals.append(val)
            # print vals
            print print_str % tuple(vals)
        sys.exit()
    def getAllSessionInfo(self):
        info("Getting Internal ID for all MR Sessions", verbosity_level = 2)
        # print str(self.session_cols)
        url = ("/data/archive/experiments?columns="+",".join(self.session_cols)+"&format=json")
        return self._sessionsFromJsonUrl(url)
    def _sessionsFromJsonUrl(self,url):
        sessions = self._getXnatJsonResponse(url)
        for session in sessions:
          if("ID" not in session):
              error("Could not get Internal ID for MR Session '" + label + "': no 'ID' in server response")
          # Make sure each "column" is populated to avoid stack traces and other
          # nasties.
          for col in self.session_cols:
              lc_col = col.lower()
              short_col = self.nons.sub('',lc_col)
              if not col in session:
                  if not lc_col in session:
                      session[lc_col] = ''
                  session[col] = session[lc_col]
              if not lc_col in session:
                  session[lc_col] = session[col]
              if not short_col in session:
                  session[short_col] = session[col]
        return sessions
    def _getXnatJsonResponse(self,url):
        (response, content) = self.geturl(url)
        obj = json.loads(content)
        if("ResultSet" not in obj):
            error("Could not get MR Session information: no 'ResultSet' in server response. Please check your input.")
        results = obj["ResultSet"]["Result"]
        if(not results):
            error("No such  MR Session found at '"+self.args.host+"'. Please make sure the MR Session is correct and that you have privileges to see it. ")
        return results
    def getSubjectInfoById(self,subject_id):
        url = "/data/archive/subjects?columns=ID,label,dob,gender,handedness,race,ethnicity,weight,height,xnat:demographicData/english,xnat:demographicData/language&ID=" + urllib2.quote(subject_id) + "&format=json"
        (response, content) = self.geturl(url)
        obj = json.loads(content)
        if("ResultSet" not in obj):
            error("Could not get info MR Session label '" + args.session_label + "': no 'ResultSet' in server response")
        results = obj["ResultSet"]["Result"]
        if(not results):
            error("Could not get info for MR Session label '" + args.session_label + "': no 'ResultSet' in server response")
        elif(len(results) > 1):
            error("Multiple Subjects returned for MR Session label '" + args.session_label + "'")
        subject = results[0];
        subject['id'] = subject['ID']
        # subject['label'] = subject['label']
        # subject['dob'] = subject['dob']
        # subject['gender'] = subject['gender']
        # subject['handedness'] = subject['handedness']
        # subject['race'] = subject['race']
        # subject['ethnicity'] = subject['ethnicity']
        # subject['weight'] = subject['weight']
        # subject['height'] = subject['height']
        if 'xnat:demographicData/english' in subject:
            subject['english'] = subject['xnat:demographicData/english']
        else:
            subject['english'] = ''
        if 'xnat:demographicData/language' in subject:
            subject['language'] = subject['xnat:demographicData/language']
        else:
            subject['language'] = ''
        return subject
    def getSeriesInfoBySession(self,session):
        url = "/data/archive/experiments/" + urllib2.quote(session['id']) + "/scans?format=json"
        (response, content) =self.geturl(url)
        obj = json.loads(content)
        if("ResultSet" not in obj):
            error("Could not get Scan information for '" + args.session_label + "': no 'ResultSet' in server response")
        series_list = obj["ResultSet"]["Result"]
        if(not series_list):
            error("No scans found for '" + args.session_label + "'")
        for series in series_list:
            series['id'] = series['ID']
            # series['quality'] = series['quality']
            # series['type'] = series['type']
            series['description'] = series['series_description']
            # series['note'] = series['note']
        if "ALL" == self.args.raw_types:
            for series in series_list:
                series['download'] = 'YES'
        else:
            series_matches = self.args.raw_types.split(',')
            for series in series_list:
                series['download'] = 'no'
                for match in series_matches:
                    if match == series['id']:
                        series['download'] = 'YES'
                    elif match == series['type']:
                        series['download'] = 'YES'
                    elif match == series['description']:
                        series['download'] = 'YES'
        session['series'] = series_list
        # print "Series 1: '"+str(session['series'])+"'"
        return series_list
    def rate(self,t0,bytes):
        mb = bytes / 1048576.0
        secs = time.time() - t0
        rate = mb/secs
        return " %0.2f MB %0.2f MB/s" % ( mb, rate )

    def getZipInFile(self,session,url):
        (request, response) = self.getRequest(url)
        count = 0
        bytes = 0
        zipFileBaseName = self.args.out_dir + "/" + session['label']
        zipFileName = zipFileBaseName + ".zip"
        while(os.path.exists(zipFileName)):
            warn("'"+zipFileName+"' already exists. Will download with different name!")
            count += 1;
            zipFileBaseName = self.args.out_dir + "/" + session['label'] + ("-copy-%03d" % (count))
            zipFileName = zipFileBaseName + ".zip"
        info("Saving to '"+zipFileName+"'")
        with open(zipFileName, 'w') as zipFile:
            try:
                #response = urllib2.urlopen(request);
                count = 0;
                while 1:
                    data = response.read(4096)
                    bytes = bytes + len(data)
                    if not data:
                        if verbosity > 0:
                            num = 50 - ((count % 5000) / 100)
                            pad = "".ljust(num," ")
                            print pad + self.rate(request.t0,bytes)+" Total"
                        break
                    zipFile.write(data);
                    count = count + 1
                    if verbosity > 0 and 0 == count % 100:
                        sys.stdout.write("#")
                        sys.stdout.flush()
                        if 0 == count % 5000:
                            print self.rate(request.t0,bytes)
                self.tock(request.t0)
            except Exception, e:
                self.tock(request.t0)
                error("Unknown Error. Can't continue because I'm not smart enough to know what to do.\nError: "+str(e))
        if self.args.legacy:
            info("Converting zip file to legacy format")
            current_zipfile = zipfile.ZipFile(zipFileName)
            legacyZipFileName = zipFileBaseName + "-legacy.zip"
            if(self.args.zip64):
              info("Using Zip64 extensions")
              legacy_zipfile = zipfile.ZipFile(legacyZipFileName, 'w', allowZip64=True)
            else:
              legacy_zipfile = zipfile.ZipFile(legacyZipFileName, 'w')
            for member in current_zipfile.namelist():
                filename = os.path.basename(member)
                source = io.BytesIO(current_zipfile.open(member).read())
                try:
                    gz_source = gzip.GzipFile(fileobj=source, mode="rb")
                    gz_source.read()
                    source = gz_source
                except IOError, e:
                    pass
                source.seek(0)
                legacy_zipfile.writestr(os.path.join(session['label'], "RAW", filename), source.read())
            legacy_zipfile.close()
            current_zipfile.close()
            os.remove(zipFileName)
            os.rename(legacyZipFileName,zipFileName)
        return (response)
    def getZipInRAMAndUnzip(self,session,url):
        (request, response) = self.getRequest(url)
        count = 0
        bytes = 0
        zipFileContents = io.BytesIO()
        try:
            response = urllib2.urlopen(request);
            count = 0;
            while 1:
                data = response.read(4096)
                bytes = bytes + len(data)
                if not data:
                    if verbosity > 0:
                        num = 50 - ((count % 5000) / 100)
                        pad = "".ljust(num," ")
                        print pad + self.rate(request.t0,bytes)+" Total"
                        break
                zipFileContents.write(data);
                count = count + 1
                if verbosity > 0 and 0 == count % 100:
                    sys.stdout.write("#")
                    sys.stdout.flush()
                    if 0 == count % 5000:
                        print self.rate(request.t0,bytes)
            self.tock(request.t0)
        except Exception, e:
            self.tock(request.t0)
            error("Unknown Error. Can't continue because I'm not smart enough to know what to do.\nError: "+str(e))
        zipFile = zipfile.ZipFile(zipFileContents)
        if self.args.new_structure is not True:
            destinationDir = os.path.join(self.args.out_dir, session['label'],'RAW')
        else:
            destinationDir = self.args.out_dir
        if not os.path.exists(destinationDir):
            os.makedirs(destinationDir)
        for member in zipFile.namelist():
            filename = os.path.basename(member)
            source = io.BytesIO(zipFile.open(member).read())
            try:
                gz_source = gzip.GzipFile(fileobj=source, mode="rb")
                gz_source.read()
                source = gz_source
            except IOError, e:
                pass
            source.seek(0)
            destination = file(os.path.join(destinationDir,filename), "wb")
            shutil.copyfileobj(source, destination)
            source.close()
            destination.close()
        return (response)
    def outputSessionDetails(self,session):
        subject = session['subject']
        series = session['series']
        print """
SUBJECT: %-30s AGE:       %-30s
GENDER:  %-30s HADEDNESS: %-30s
RACE:    %-30s ETHNICITY: %-30s
ENGLISH: %-30s LANGUAGE:  %-30s
WEIGHT:  %-30s HEIGHT:    %-30s
SESSION: %s (Accession ID %s)
""" % (subject['label'] + " ("+subject['id']+")", self.getage(subject['dob'],session['date']),
       subject['gender'], subject['handedness'],
       subject['race'], subject['ethnicity'],
       subject['english'], subject['language'],
       subject['weight'], subject['height'],
       session['label'], session['id'])
        print "%6s %4s %-8s %-12s %-32s %s" % ('downld','num','quality','type','series','note');
        for series in session['series']:
            # print "Series: '" +str(series)+"'"
            # print "series['download']='"+str(series['download'])+"'"
            # print "series['ID']='"+str(series['ID'])+"'"
            # print "series['quality']='"+str(series['quality'])+"'"
            # print "series['type']='"+str(series['type'])+"'"
            # print "series['description']='"+str(series['description'])+"'"
            # print "series['note']='"+str(series['note'])+"'"
            print "%6s %4s %-8s %-12s %-32s %s" % (series['download'],series['ID'],series['quality'],series['type'], series['description'],series['note'])
        print ""
    def downloadSeriesToZipFile(self,session):
        url = "/data/archive/experiments/" + urllib2.quote(session['id']) + "/scans/" + urllib2.quote(self.args.raw_types) + "/files?format=zip"
        return self.getZipInFile(session,url)
    def downloadSeriesToDir(self,session):
        for series in session['series']:
            if 'YES' == series['download']:
                info("Scan %4s %-12s %-24s from %s" %
                     (series['id'],series['type'],series['description'], "'"+session['label']+"'"))
                url = ("/data/archive/experiments/" + urllib2.quote(session['id']) + "/scans/" +
                       series['id'] + "/files?format=zip")
                (response) = self.getZipInRAMAndUnzip(session,url)
        info("Done with '"+session['label']+"'.")
def main():
  ## --- arguments
  parser = argparse.ArgumentParser(description="ArcGet: retrieve imaging data from XNAT")
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument("--host", "-host", required=False, help="XNAT URL")
  group.add_argument("-a", "--alias", required=False, help="XNAT config file alias")
  parser.add_argument("-l", "--legacy", required=False, action="store_true", help="Return legacy XNAT 1.4 zipfile and directory structure")
  parser.add_argument("--new-structure", required=False, default=False,action="store_true", help="Don't create SessionID/RAW directory, output files in current working directory or --out-dir")
  # parser.add_argument("--clean-names", required=False, action="store_true", help="Make all file have the SESSIONID_series_SERIESNUM_file_FILENUM.dcm format")
  parser.add_argument("-u", "--username", required=False, help="XNAT username")
  parser.add_argument("-p", "--password", required=False, help="XNAT password")
  parser.add_argument("-s", "--session-label", action="append", required=False,dest="session_labels",help="MR Session label")
  parser.add_argument("-r", "--raw-types", required=False, default="ALL", help="raw scan types or numbers e.g. 1,MEMPRAGE,21,22,DSI")
  parser.add_argument("-o", "--out-dir", default=".", required=False, help="output directory")
  parser.add_argument('-q', "--quiet", nargs=0, action=ArgParseSubAction, dest='quiet', help="Decrease verbosity by 1. Can be used several times.")
  parser.add_argument('-v', "--verbose", nargs=0, action=ArgParseAddAction, dest='verbose', help="Increase verbosity by 1. Can be used several times.")
  parser.add_argument('--zip64', required=False, action="store_true", help="Use Zip64 extensions when creating zip archives")
  parser.add_argument("-W", "--no-warn", required=False, action="store_true", help="Don't show me the annoying warning about downloading everything again. I like wasting bandwidth.")
  parser.add_argument("--show-all",  action='store_true', dest='show_all', help="Show some information on all available sessions and exit.")
  (args,sessions) = parser.parse_known_args()
  if(args.password and not args.legacy):
    error("DO NOT put passwords on the command line unless absolutely necessary!! --password only allowed with --legacy")
  # print "================================================================="
  # print "Before:"
  # for arg in args:
  #     print arg + "='"+str(args[arg])+"'"
  # print "================================================================="
  ## --- read username and password from XNAT config file
  config_file = os.path.expanduser("~/.xnat_auth")
  if(not os.path.isfile(config_file)):
    info("No config file found: " + config_file)
    if(args.alias):
      error("You cannot specify an --alias without a config file")
  else:
    info("Reading config file: " + config_file, verbosity_level = 3)
    xml = etree.parse(os.path.expanduser(config_file))
    if(args.alias):
      if(not args.alias.isalnum()):
        error("--alias must be alphanumeric", parser=parser)
      ## --- get host
      args.host = xml.xpath("/xnat/" + args.alias + "/url/text()")
      if(args.host):
        args.host = args.host[0]
      ## --- get username
      if(not args.username):
        args.username = xml.xpath("/xnat/" + args.alias + "/username/text()")
        if(args.username):
          args.username = args.username[0]
      ## --- get password
      args.password = xml.xpath("/xnat/" + args.alias + "/password/text()")
      if(args.password):
        args.password = args.password[0]
    elif(args.host):
      ## --- get username
      if(not args.username):
        args.username = xml.xpath("/xnat/*[url='" + args.host + "']/username/text()")
        if(args.username):
          args.username = args.username[0]
      ## --- get password
      if(not args.password):
        args.password = xml.xpath("/xnat/*[url='" + args.host + "']/password/text()")
      	if(args.password):
        	args.password = args.password[0]
  ## --- prompt for host, username, password if necessary
  if(sys.stdin.isatty()):
    if(not args.host):
      args.host = raw_input("Enter host: ")
    if(not args.username ):
      args.username = raw_input("Enter username: ")
    if(not args.password):
      args.password = getpass.getpass("Enter password: ")
  ## --- final check
  if(not args.host):
    error("Could not retrieve a host from config file or command line")
  if(not args.username):
    error("Could not retrieve a username from config file or command line")
  if(not args.password):
    error("Could not retrieve a password from config file or command line")
  ## --- strip any slashes from right side of host
  args.host = str(args.host).rstrip("/")
  if not os.path.exists(args.out_dir):
      os.makedirs(args.out_dir)
  args.out_dir = os.path.abspath(args.out_dir)
  info("Saving output to '"+args.out_dir+"'")
  if args.session_labels is None:
      args.session_labels = []
  if len(sessions) > 0:
      args.session_labels.extend(sessions)
  # print str(args)
  # sys.exit()
  # log(str(args))

  # print "================================================================="
  # print "After:"
  # for arg in args:
  #     print arg + "='"+str(args[arg])+"'"
  # print "================================================================="
  # sys.exit()
  arcget = ArcGet(args)
  if args.show_all:
      arcget.getShowAllSessionInfo()
      sys.exit()
  sessions = arcget.getSessionInfoByLabels(args.session_labels)

  for session in sessions:
      subject = arcget.getSubjectInfoById(session['subject_id'])
      if "ALL" == args.raw_types:
          if not args.no_warn:
              warn("+-------------------------------------------------------------------+",verbosity_level=-1)
              warn("|               --==>> WARNING: READ CAREFULLY <<==--               |",verbosity_level=-1)
              warn("+-------------------------------------------------------------------+",verbosity_level=-1)
              warn("| By not specifying which scans/series to download from the session |",verbosity_level=-1)
              warn("| you will be downloading EVERYTHING, including report files, text  |",verbosity_level=-1)
              warn("| files, pictures, and EVERY SINGLE scan. If you don't REALLY NEED  |",verbosity_level=-1)
              warn("| it all you are saying that you really do want to waste EVERYONE's |",verbosity_level=-1)
              warn("| space, processing power, time, and slow down XNAT, the cluster    |",verbosity_level=-1)
              warn("| etc, etc. So DON'T DO IT. Use the -r or --raw-types option, for   |",verbosity_level=-1)
              warn("| example, to get the first, third and all BOLD scans, use:         |",verbosity_level=-1)
              warn("|    --raw-types 1,3,BOLD                                           |",verbosity_level=-1)
              warn("+-------------------------------------------------------------------+",verbosity_level=-1)
      session['subject'] = subject
      series_list  = arcget.getSeriesInfoBySession(session)
      if verbosity > 0:
          arcget.outputSessionDetails(session)
      if(args.legacy):
          info("Creating legacy XNAT zipfile")
          arcget.downloadSeriesToZipFile(session)
      else:
          info("Getting Data...")
          arcget.downloadSeriesToDir(session)

def error_noexit(message=""):
  '''
  Error message
  @param message: str
  @param status: int
  '''
  global progname
  print(progname+": ERROR: " + str(message))
def error(message="", status=1, parser=None):
  '''
  Error message
  @param message: str
  @param status: int
  '''
  global verbosity
  global progname
  error_noexit(message)
  if(parser):
    parser.print_help()
  if(verbosity > 2):
    traceback.print_exc(file=sys.stdout)
  sys.exit(status)
def info(message="",verbosity_level=1):
  '''
  Information message
  @param message: str
  '''
  global verbosity
  if(verbosity >= verbosity_level and message):
      print(str(message))
warnings=[]
def warn(message="",verbosity_level=0):
  '''
  Warning message
  @param message: str
  '''
  global verbosity,warnings
  global progname
  warnings.append(message)
  if(verbosity >= verbosity_level):
      print(progname+": WARNING: " + str(message))

def repeatWarnings():
  '''
  Repeat warning for good measure when the script exits
  '''
  global verbosity,warnings
  global progname
  if(verbosity >= 0 and len(warnings) > 0):
      print progname + ": ================================================================="
      print progname + ": SUMMARY of WARNINGS:"
      for warning in warnings:
          print(progname+": WARNING: " + str(warning))
# Do after defining repeatWarnings()
import atexit
atexit.register(repeatWarnings)
        
def debug(message="", verbosity_level = 4):
  '''
  Debug message
  @param message: str
  '''
  global verbosity
  global progname
  if(verbosity >= verbosity_level):
      print(progname+": DEBUG: " + str(message))
class ArgParseAddAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        global verbosity
        #print "values='"+str(values)+"'"
        # print 'values: {v!r}'.format(v=values)
        if values==None:
            values='1'
        elif isinstance(values, list):
            values='1'
        try:
            values=int(values)
        except ValueError:
            values=values.count('v')+1
        verbosity = verbosity + values
        setattr(args, self.dest, values)
class ArgParseSubAction(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        global verbosity
        print "values='"+str(values)+"'"
        # print 'values: {v!r}'.format(v=values)
        if values==None:
            values='1'
        elif isinstance(values, list):
            values='1'
        try:
            values=int(values)
        except ValueError:
            values=values.count('v')-1
        verbosity = verbosity - values
        setattr(args, self.dest, values)
class ArgParseAddActionOld(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        global verbosity
        # print 'values: {v!r}'.format(v=values)
        if values==None:
            values='1'
        try:
            values=int(values)
        except ValueError:
            values=values.count('v')+1
        verbosity = verbosity + values
        setattr(args, self.dest, values)
class ArgParseSubActionOld(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        global verbosity
        # print 'values: {v!r}'.format(v=values)
        if values==None:
            values='1'
        try:
            values=int(values)
        except ValueError:
            values=values.count('v')+1
        verbosity = verbosity - values
        setattr(args, self.dest, values)

if __name__ == "__main__":
    main()
