# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os
import uuid
import logging
import jinja2
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace as Workspace

class BaseModule:
    def __init__(self,config,version,name):
        self.config = config
        if "SDK_CALLBACK_URL" in os.environ:
            self.callback_url = os.environ['SDK_CALLBACK_URL']
            self.dfu = DataFileUtil(self.callback_url)
        self.version = version
        self.name = name
        self.scratch_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
            level=logging.INFO)
        self.clear_context()
        self.report_html = None
    
    def validate_args(self,params,required,defaults):
        for item in required:
            if item not in params:
                raise ValueError('Required argument '+item+' is missing!')
        for key in defaults:
            if key not in params:
                params[key] = defaults[key]
        return params
    
    def clear_context(self):
        self.report_info = None
        self.ctx = None
        self.output_type = None
        self.output_id = None
        self.wsclient = None
    
    def finalize_call(self,output):
        if self.report_info != None:
            output['report_name'] = self.report_info['name']
            output['report_ref'] = self.report_info['ref']
        if self.workspace != None:   
            output['workspace_name'] = self.workspace
            output['ws'] = self.workspace
        if self.output_type != None:
            output['type'] = self.output_type
            output['obj'] = self.output_id
        return output

    def initialize_call(self,ctx,workspace=None,output_type = None,output_id = None):
        self.clear_context()
        self.workspace = workspace
        self.ctx = ctx
        self.output_type = output_type
        self.output_id = output_id
        self.objects_created = []
        self.wsclient = Workspace(self.config["workspace-url"], token=self.ctx['token'])
    
    def add_created_object(self,ref,description):
        self.objects_created.append({"ref":ref,"description":description})

    def create_report(self,context,template_file=None,height=500):
        html_report_folder = os.path.join(self.scratch_folder, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        
        with open(os.path.join(html_report_folder, 'view.html'), 'w') as f:
            self.report_html = self.build_report(context,template_file)
            f.write(self.report_html)

        report_shock_id = ""
        if self.config["save_report_to_kbase"] == "1":
            report_shock_id = self.dfu.file_to_shock({'file_path': html_report_folder,'pack': 'zip'})['shock_id']

        html_output = {
            'name' : 'view.html',
            'shock_id': report_shock_id
        }
        report_params = {
            'objects_created': self.objects_created,
            'workspace_name': self.workspace,
            'html_links': [html_output],
            'direct_html_link_index': 0,
            'html_window_height': height,
            'report_object_name': self.name + '_report_' + str(uuid.uuid4())
        }
        if self.config["save_report_to_kbase"] == "1":
            report = KBaseReport(self.callback_url, token=self.ctx['token'])
            self.report_info = report.create_extended_report(report_params)
        return self.report_html
        
    def build_report(self,context,template_file=None):
        if template_file == None:
            template_file = self.config["template_file"]
        # Directory this file is in
        array = template_file.split("/")
        filename = array.pop()
        template_dir = "/".join(array)
        env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(template_dir),
                autoescape=jinja2.select_autoescape(['html', 'xml']))
        # Return string of html
        return env.get_template(filename).render(context)