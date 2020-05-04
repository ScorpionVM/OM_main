import json, math

from components import PRINT, INPUT
PRT = PRINT()
INP = INPUT()

# tools
def sprijin(moment):
    if moment == "Mta":
        return "pe placuta din bronz cu suprafata marita"
    elif moment == "Mtb":
        return "pe placuta din otel de suprafata sferica"
    else:
        return "otel cu suprafata sferica cu raza marita"

def json_read(path):
    with open(path,"r") as f:
        data = json.load(f)
        f.close()
    return data

def write_json(path, data):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
        f.close()

def absolute(x):
    if x < 0:
        return x * (-1)
    else:
        return x

# parser de eroare la extragere
def definire_date(etapa, path, diam):
    PRT.pcol("[-] Eroare de extragere a datelor pentru diametru [ ø %.2f ] " % diam)
    PRT.header(":: Introduceti datele in baza de date pentru un diametru nou ::")
    
    add_new = {}

    if etapa == 1:
        PRT.pcol("[!] Adaugati la ajustaje pentru diametrul ø %.2f :" % diam)
        limite = INP.icol("[?] Introduceti limitele de diametru x-y : ")
        add_new[limite] = {}
        while True:
            if 'n' == INP.icol("[?] Add new value [Alezaj] = [y/N] -> "): 
                break
            
            abat_alz = INP.icol("[?] Abatere alezaj -> ")
            add_new[limite][abat_alz] = {
                "main": {
                    "AS" : int(INP.icol("[?] %s AS -> " % abat_alz)),
                    "AI" : int(INP.icol("[?] %s AI -> " % abat_alz))
                }, 
                "date": {} 
            }
            PRT.pcol("[!] Pentru alezaj %s / adaugati arbori :" % abat_alz)
            while True:
                if 'n' == INP.icol("[?] Add new value [Arbori] = [y/N] -> "): 
                    break

                arbori = INP.icol("[?] %s / arbore -> " % abat_alz)
                add_new[limite][abat_alz]["date"][arbori] = {
                    "as": int(INP.icol("[?] %s as -> " % arbori)),
                    "ai": int(INP.icol("[?] %s ai -> " % arbori))
                }

            if 'y' == INP.icol("[?] Add extended [y/N] -> "):
                add_new[limite][abat_alz]["date"]["extended"] = {}
                while True:
                    if 'n' == INP.icol("[?] Add new value [Extended]= [y/N] -> "): 
                        break
                    
                    extlimite = INP.icol("[?] Introduceti limitele extinse de diametru x-y : ")
                    add_new[limite][abat_alz]["date"]["extended"][extlimite] = {}
                    while True:
                        if 'n' == INP.icol("[?] Add new value [Extended/Arbori] [y/N] -> "): 
                            break
                        
                        arbori = INP.icol("[?] %s / arbore -> " % abat_alz)
                        add_new[limite][abat_alz]["date"]["extended"][extlimite][arbori] = {
                            "as": int(INP.icol("[?] %s as -> " % arbori)),
                            "ai": int(INP.icol("[?] %s ai -> " % arbori))
                        }
            else:
                add_new[limite][abat_alz]["date"]["extended"] = None
        
        Block = json_read("./AbateriAjustaje.json")
        Block.update(add_new)
        write_json("./AbateriAjustaje.json", Block)


            
    elif etapa == 2.1:
        PRT.pcol("[!] Adaugati la dimensiunile profil patratic pentru diametrul ø %.2f :" % diam)
        limite = INP.icol("[?] Introduceti limitele de diametru x-y : ")

    elif etapa == 2.2:
        PRT.pcol("[!] Adaugati la dimensiunile profil trapezoidal pentru diametrul ø %.2f :" % diam)
        limite = INP.icol("[?] Introduceti limitele de diametru x-y : ")

# Etapa 1
def getDiametru(diam):
    json_path = "./AbateriAjustaje.json"
    try:
        Diametru = json_read(json_path)
        for key in Diametru:
            limits = key.split('-')
            if float(limits[0]) < diam and diam < float(limits[-1]):
                diams = Diametru[key]

        extended = None

        for key in diams:
            if diams[key]["date"]["extended"] is not None and extended is None:
                for ext_key in diams[key]["date"]["extended"]:
                    limits = ext_key.split('-')
                    if float(limits[0]) < diam and diam < float(limits[-1]):
                        diams[key]["date"].update(diams[key]["date"]["extended"][ext_key])

            del diams[key]["date"]["extended"]

        return diams
    except Exception:
        #definire_date(1,json_path, diam)
        print(Exception)

class OteluriLaminateCarbon:
    Marca = {
        "OLC 15": {"Rp02": 350},
        "OLC 20": {"Rp02": 250},
        "OLC 45": {"Rp02": 360},
        "OLC 60": {"Rp02": 400},
    }

# Etapa 2
def getDiametru_d1_FiletPatrat(diam):
    json_path = "./DimensiuneProfilFiletPatrat.json"
    # try:
    Diametru = json_read(json_path)
    # per diametru
    for key in Diametru:
        limits = key.split("-")
        if float(limits[0]) < diam and diam <= float(limits[-1]):
            block = Diametru[key]
            
    for i in range(len(block["intv"])-1):
        if block["intv"][i] <= diam and diam <= block["intv"][i+1]:
            d1 = block["intv"][i+1]
            diametre = {
                "d": d1 + block["p"],
                "d1": d1,
                "d2": d1 + block["p"] * 0.5,
                "D2": d1 + block["p"] * 0.5,
                "p": block["p"]
            }
    return diametre
    #except Exception:
    #    definire_date(2.1, json_path, diam)

def getDiametru_d1_FiletTrapez(diam):
    json_path = "./DimensiuneProfilFiletTrapezoidal.json"
    #try:
    Diametru = json_read(json_path)

    for key in Diametru:
        limits = key.split("-")
        if float(limits[0]) < diam and diam <= float(limits[-1]):
            block = Diametru[key]
    
    for i in range(len(block["intv"])-1):
        if block["intv"][i] < diam and diam <= block["intv"][i+1]:
            d1 = block["intv"][i+1]
            d = math.trunc(d1 + block["p"] + 1)
            diametre = {
                "d": d,
                "d1": d1,
                "d2": d - block["p"] * 0.5,
                "D2": d - block["p"] * 0.5,
                "p": block["p"]
            }
        
    return diametre
    #except Exception:
    #    definire_date(2.2, json_path, diam)


def getDiametru_di(diam):
    recomandari = {
        "6-16": [2.12, 2.65],
        "17-32": [3.35, 4.25],
        "33-56": [5.30, 6.70],
        "57-80": [8.50, 10.6],
        "81-120": [16, 18]
    }

    for key in recomandari:
        limits = key.split("-")
        if float(limits[0]) <= diam and diam <= float(limits[-1]):
            int_d = float(limits[-1]) - float(limits[0])
            int_di = recomandari[key][-1] - recomandari[key][0]
            cf = int_di/int_d
            di = recomandari[key][0] + cf

    return di
    