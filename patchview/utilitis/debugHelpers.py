# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 14:58:25 2021

@author: MHu
"""

from inspect import getframeinfo, stack
import pdb

def debugInfo(message, setBreak=True):
    caller = getframeinfo(stack()[1][0])
    print(f"\n{caller.filename}: line {caller.lineno}, \n\t{message}")
    if setBreak:
        pdb.set_trace()
