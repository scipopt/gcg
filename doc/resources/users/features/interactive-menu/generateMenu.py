#!/usr/bin/env python3
import subprocess
import signal
import os


def getItem(command = "help", level = 0):
    quitCommand = "-c 'quit'"
    goBackCommand = ""
    execCommand = ""
    for com in command.split():
        execCommand += "-c '" + com + "'"
    print("Level:" + str(level))
    i = 0
    while i < level:
        goBackCommand += '-c ".."'
        i += 1

    #print(f'Executing "{command}"')
    proc = subprocess.Popen([f'./../../../../../bin/gcg {execCommand} {goBackCommand} {quitCommand} | grep -A100 -m1 "user parameter file" | tail -n+4 | sed "s/^  //g" | sed "/\\n/d"'], shell=True, stdout=subprocess.PIPE,universal_newlines=True)

    try:
        outs, errs = proc.communicate(timeout=15)
    except TimeoutExpired:
        proc.kill()
        outs, errs = proc.communicate()

    outs = str(outs).split("\n")
    return outs

def getSubmenu(menu, level = 1):
    menu_temp = []
    for i in range(len(menu)):
        menu_temp.append(menu[i])
        if menu[i].startswith("<"):
            to_insert = getItem(menu[i].split('<')[1].split('>')[0],level=level)
            for item in to_insert:
                print("Inserting: " + item)
                menu_temp.append(item)
                menu.append(getSubmenu(item, level = level+1))
        i += 1
    return menu_temp

def main():
    # get help first
    menu = getItem(command = "help", level = 0)
    menu = getSubmenu(menu, level = 1)
    print(menu)


if __name__ == '__main__':
    main()
