#!/usr/bin/env python

import sys, os
import argparse
import re
import os.path
import subprocess
import errno

def mkdir_p(path):
    try:
        mkdir_command=" ".join(["mkdir -p", path])
        os.system(mkdir_command)

    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

current_dir=os.getcwd()
main_scripts_directory=os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('-gwas',help='Complete path with GWAS file',action='store')
parser.add_argument('-priority_regions',help='Complete path with Priority Regions file',action='store')
parser.add_argument('-out_dir',help='OPTIONAL:Complete path for output files',action='store')
parser.add_argument('-prefix',help='OPTIONAL:Prefix for output files',action='store')
parser.add_argument('-task',help='Choose to run: sFDR or wFDR analysis',action='store' )
parser.add_argument('-column', help='OPTIONAL: Weight column in Priority Regions for transformation',type=int)
args=parser.parse_args()
error=1;

if args.gwas and not args.priority_regions:
    parser.error("Please add Priority Regions file")
elif not args.gwas and args.priority_regions:
    parser.error ("Please add GWAS file")
elif not args.gwas and not args.priority_regions:
    parser.error ("Please add GWAS and Priority Regions file")
else:
    error=0;
    print ("INPUT FILES: GWAS and Priority Regions")
if error==0:
    if args.out_dir:
        if args.prefix:
            print ("Output directory and file name introduced correctly")
else:
            error=1
            sys.exit("Please add missing arguments")

if args.gwas:
    gwas_dir=os.path.dirname(args.gwas)
    check_dir=os.path.exists(gwas_dir)
    if check_dir==False:
        sys.exit("INCORRECT GWAS DIRECTORY")
    else:
        print("GWAS DIRECTORY: "),gwas_dir;
    os.chdir(gwas_dir)
    gwas_file=os.path.basename(args.gwas)
    check_file=os.path.isfile(gwas_file)
    if check_file==False:
        sys.exit("INCORRECT GWAS FILENAME")
    else:
        print("GWAS FILENAME: "),gwas_file;
    gwas_a_error=0
    gwas_b_error=0
    gwas_c_error=0
    gwas_d_error=0
    gwas_errors=1
    x=0
    with open(gwas_file) as file:
        for line in file:
            x=x+1;
            if x==6:
                break
            else:
                fields = line.split( "\t" )
                gwas_bed = (fields[0], int(fields[1]), int(fields[2]), fields[3] )
                column_a=(fields[0])
                column_b=(fields[1])
                column_c=(fields[2])
                column_d=(fields[3])
                r = re.compile("([a-zA-Z]+)([0-9]+)")
                m = r.match(column_a)
                n = r.match(column_d)
                if m.group(1) == 'chr':
                    gwas_a_error=gwas_a_error
                else:
                    gwas_a_error=gwas_a_error+1
                if column_d.isalnum() == True:
                    gwas_d_error=gwas_d_error
                else:
                    gwas_d_error=gwas_d_error+1
                
                column_a=m.group(2)
                column_d=n.group(2)
                a_ask=column_a.isdigit()
                b_ask=column_b.isdigit()
                c_ask=column_c.isdigit()
                d_ask=column_d.isdigit()
                if a_ask==True:
                    gwas_a_error=gwas_a_error
                else:
                    gwas_a_error=gwas_a_error+1
                if b_ask==True:
                    gwas_b_error=gwas_b_error
                else:
                    gwas_b_error=gwas_b_error+1
                if c_ask==True:
                    gwas_c_error=gwas_c_error
                else:
                    gwas_c_error=gwas_c_error+1
                if d_ask==True:
                    gwas_d_error=gwas_d_error
                else:
                    gwas_d_error=gwas_d_error+1
    if gwas_a_error!=0:
        print("Warning: Your GWAS file doesn't have a chr# written in the first column")
    if gwas_d_error!=0:
        print("Warning: Your GWAS file doesn't have a rs# written in the fourth column")
    if gwas_a_error==0 and gwas_b_error==0 and gwas_c_error==0 and gwas_d_error==0:
        gwas_errors=0;
        print "Gwas file: CORRECT FORMAT"
    else:
        print "INCORRECT FORMAT"
        sys.exit ("GWAS files is not in the correct format")

if args.priority_regions:
    p_regions_dir=os.path.dirname(args.priority_regions)
    check_dir=os.path.exists(p_regions_dir)
    if check_dir==False:
        sys.exit("INCORRECT PRIORITY REGIONS DIRECTORY")
    else:
        print("PRIORITY REGIONS DIRECTORY: "),p_regions_dir;
    os.chdir(p_regions_dir)
    p_regions_file=os.path.basename(args.priority_regions)
    check_file=os.path.isfile(p_regions_file)
    if check_file==False:
        sys.exit("INCORRECT PRIORITY REGIONS FILENAME")
    else:
        print("PRIORITY REGIONS FILENAME: "),p_regions_file;
    regions_a_error=0
    regions_b_error=0
    regions_c_error=0
    regions_errors=1
    y=0
    with open(p_regions_file) as file:
        for line in file:
            y=y+1;
            if y==6:
                break
            else:
                fields = line.split( "\t" )
                regions_bed = (fields[0], int(fields[1]), int(fields[2]) )
                column_a=fields[0]
                column_b=(fields[1])
                column_c=(fields[2])
                r = re.compile("([a-zA-Z]+)([0-9]+)")
                m = r.match(column_a)
                if m.group(1) == 'chr':
                    regions_a_error=regions_a_error
                else:
                    regions_a_error=regions_a_error+1
                column_a=m.group(2)
                a_ask=column_a.isdigit()
                b_ask=column_b.isdigit()
                c_ask=column_c.isdigit()
                if a_ask==True:
                    regions_a_error=regions_a_error
                else:
                    regions_a_error=regions_a_error+1
                if b_ask==True:
                    regions_b_error=regions_b_error
                else:
                    regions_b_error=regions_b_error+1
                if c_ask==True:
                    regions_c_error=regions_c_error
                else:
                    regions_c_error=regions_c_error+1
    if regions_a_error==0 and regions_b_error==0 and regions_c_error==0:
        regions_errors=0;
        print "Priority Regions file:CORRECT FORMAT"
    else:
        print "INCORRECT FORMAT"
        sys.exit ("Priority regions file is not in the correct format")

os.chdir(current_dir)

if args.out_dir:
    out_dir=args.out_dir
    mkdir_p(out_dir)
    full_out_dir=os.path.join(current_dir,out_dir)
    print "OUTPUT DIRECTORY: ", full_out_dir;
if not args.out_dir:
    out_dir='Results_Directory'
    mkdir_p(out_dir)
    full_out_dir=os.path.join(current_dir,out_dir)
    print "OUTPUT DIRECTORY created: ", full_out_dir;

if args.prefix:
    prefix=args.prefix
    prefix_path=os.path.join(full_out_dir,prefix)
    ask=os.path.exists (prefix_path)
    if ask==False:
        print "PREFIX NAME: ", prefix;
    
    else:
        print "PREFIX NAME: ", prefix;
        sys.exit("WARNING: The name of the file is already asigned. This option will overwrite your files. If you wish to do so, enter the same name")
if not args.prefix:
    prefix="result"
    prefix_path=os.path.join(full_out_dir,prefix)
    ask=os.path.exists (prefix_path)
    if ask==False:
        print "PREFIX NAME created: ", prefix;
    else:
        print "PREFIX NAME created: ", prefix;
        sys.exit("WARNING: The name of the file is already asigned. This option will overwrite your files. If you wish to do so, enter the same name")

if args.task=='sFDR':
    r_script='run.hd.gwas.singlesnp.sfdr.R'
    col="NA"
    dim="NA"

if args.task=='wFDR':
    r_script='run.hd.gwas.singlesnp.wfdr.R'
    col=args.column
    col=int(col)
    print("Defining dimension of GWAS and Priority Regions files...")
    os.chdir(gwas_dir)
    with open(gwas_file) as file:
        for line in file:
            fields = line.split( "\t" )
            dim=len(fields)
    os.chdir(p_regions_dir)
    with open(p_regions_file) as file:
        for line in file:
            pfields = line.split( "\t" )
        dim_priority=len(pfields)
        dim_priority=int(dim_priority)
    if col > dim_priority:
        sys.exit("This input exceeds the size of the Priority Regions file. Enter all the information once again")
    col=str(col)
    dim=str(dim)
    os.chdir(current_dir)
    if not args.column:
        sys.exit("There is no input for the weight column. Please add missing argument")

if not args.task:
        sys.exit("Task is not specified. Please repeat complete command")


################
if gwas_errors==0 and regions_errors==0 and error==0:
    rdirectory=main_scripts_directory
    os.chdir(current_dir)
    path_r_gwas=os.path.join(rdirectory, r_script)
    print "Path to R script: ",path_r_gwas;
    answer=os.path.isfile(path_r_gwas)
    if answer==True:
        rcommand=['Rscript',path_r_gwas ,args.gwas,args.priority_regions,out_dir,prefix, main_scripts_directory, args.task, col,dim]

        rcommand_string=" ".join(rcommand)
        print (rcommand_string)
        process = subprocess.Popen(rcommand)

    else:
        print "ERROR when running HD GWAS instructions"
        sys.exit ("Rscript for the analysis is not available ")

