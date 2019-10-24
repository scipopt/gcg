#!/usr/bin/env python3
import subprocess
import signal
import os


def getSubmenu(command = "help", level = 0):
    print("Command to be executed: " + command)
    quitCommand     = " -c 'quit'"
    goBackCommand   = " -c '"
    execCommand     = " -c '"
    for com in command.split():
        execCommand += com + " "
    execCommand += "'"

    print("Level:" + str(level))
    i = 0
    while i < level:
        goBackCommand += " .."
        i += 1
    goBackCommand += "'"

    #print(f'Executing "{command}"')
    execstring = "./../../../../../bin/gcg {} {} {} | grep -A100 -m1 'user parameter file' | tail -n+4 | sed 's/^  //g' | sed '/\\n/d'".format(execCommand, goBackCommand, quitCommand)
    print(execstring)
    proc = subprocess.Popen([str(execstring)], shell=True, stdout=subprocess.PIPE,universal_newlines=True)
    try:
        outs, errs = proc.communicate(timeout=15)
    except TimeoutExpired:
        proc.kill()
        outs, errs = proc.communicate()

    outs = str(outs).split("\n")
    return outs

def getMenu(menu, level = 1, previousCmd = ""):
    print("Trying to get submenu: " + str(menu))
    menu_temp = []
    # iterate through current menu
    for i in range(len(menu)):
        print("Checking item " + menu[i])
        if menu[i].startswith("<no options available>"):
            print("========================================")
            print("= WARNING: No options for menu {} =".format(menu[i]))
            print("========================================")
            continue
        elif menu[i].startswith("<set>") and level == 1:
            print("========================================")
            print("= INFORMATION: Skipping <set> submenu. =")
            print("========================================")
            continue
        elif menu[i] == '':
            continue
        else:
            menu_temp.append(menu[i])
        # check whether there are items that have a submenu
        if menu[i].startswith("<"):
            for item in getMenu(getSubmenu(previousCmd + menu[i].split('<')[1].split('>')[0],level=level), level = level+1, previousCmd = previousCmd + menu_temp[-1].split('<')[1].split('>')[0]+ " "):
                menu_temp.append(item)
        i += 1
    return menu_temp

def main():
    # get help first
    menu = getSubmenu(command = "help", level = 0)
    menu = getMenu(menu, level = 1)
    for entry in menu:
        print(entry + "\n")


if __name__ == '__main__':
    main()
