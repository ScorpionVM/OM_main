import re
from colorama import init, Fore, Style

init()

class PRINT:
    def pcol(self, data):
        data = data.split(' ')
        symb = data[0]
        temp = " ".join(data[1:])
        value = " " + temp.replace('[', Fore.LIGHTCYAN_EX+'['+Style.RESET_ALL).replace(']', Fore.LIGHTCYAN_EX+']'+Style.RESET_ALL)

        if symb == '[!]':
            print('\n' + Fore.LIGHTYELLOW_EX + symb + Style.RESET_ALL + value)
        elif symb == '[+]':
            print(Fore.LIGHTGREEN_EX + symb + Style.RESET_ALL + value)
        elif symb == '[*]':
            print(Fore.LIGHTBLUE_EX + symb + Style.RESET_ALL + value)
        elif symb == '[-]':
            print(Fore.LIGHTRED_EX + symb + Style.RESET_ALL + value)
        elif symb == '-':
            print(Fore.LIGHTCYAN_EX + symb + Style.RESET_ALL + value)
        elif symb[-1] == ')':
            print("\n"+Fore.LIGHTCYAN_EX + " "+ symb[:-1] + Style.RESET_ALL + Fore.LIGHTRED_EX + ')' + Style.RESET_ALL + value+"\n")
        else:
            print(value)


    def header(self, data):
        data = data.replace('[', Fore.LIGHTCYAN_EX+'['+Style.RESET_ALL).replace(']', Fore.LIGHTCYAN_EX+']'+Style.RESET_ALL).replace('::', Fore.LIGHTCYAN_EX+'::'+Style.RESET_ALL)
        print('\n'+data+'\n')

class INPUT:
    def icol(self, data):
        symb = Fore.LIGHTCYAN_EX+"[?] "+Style.RESET_ALL
        data = data.replace('[', Fore.LIGHTCYAN_EX+'['+Style.RESET_ALL).replace(']', Fore.LIGHTCYAN_EX+']'+Style.RESET_ALL).replace('=', Fore.LIGHTCYAN_EX+'='+Style.RESET_ALL)
        x = input(symb+data)
        return x