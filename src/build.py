from fabricate import *

sources = ['test', 'eval']

cpp_flags = '-Wall -std=c++11'.split()
cpp_flags.append('-I./eigen')
link_flags = []

def build():
    compile()
    link()

def compile():
    for source in sources:
        run('g++', '-c', source+'.cpp', cpp_flags)

def link():
    objects = [s+'.o' for s in sources]
    run('g++', '-o', 'test', objects)

def clean():
    autoclean()

main()
