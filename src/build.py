from fabricate import *

sources = ['test']

def build():
    compile()
    link()

def compile():
    for source in sources:
        run('g++', '-c', source+'.cpp')

def link():
    objects = [s+'.o' for s in sources]
    run('g++', '-o', 'test', objects)

def clean():
    autoclean()

main()
