#!/bin/bash
 grep "# " README.md | sed 's/###/\-\ [\ ]/g' | sed '/Module/s/$/ Checklist/' | sed '/Prerequisite/d' | sed 's/Prof. David Minton/[YOUR NAME]/'> module2-checklist.md
