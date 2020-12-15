#!/bin/bash
jupyter nbconvert\
   --to html\
   --output-dir docs\
   --TemplateExporter.exclude_input=True\
   --TemplateExporter.exclude_output_prompt=True\
   $1
