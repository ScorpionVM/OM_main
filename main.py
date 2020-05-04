#!/usr/bin/python3

# from py lib
import math
from tabulate import tabulate
from colorama import Fore, Style

# from custom path
from date_stabilite import getDiametru, getDiametru_d1_FiletPatrat, getDiametru_d1_FiletTrapez, getDiametru_di, OteluriLaminateCarbon as OLC, sprijin, absolute
from components import PRINT, INPUT
PRT = PRINT()
INP = INPUT()

# start main block
if 'n' == INP.icol("Start test mode [y/N] : "):
    PRT.header(":: Datele initiale ::")
    d = float(INP.icol("d [mm] = "))
    d1 = 0
    d2 = float(INP.icol("d2 [mm] = "))
    Mt = float(INP.icol("Mt [N*mm] = "))
    l = float(INP.icol("l [mm] = "))
else:
    d=53; d1=0; d2=89; Mt=115000; l=43
    PRT.pcol("[+] d = %d [mm]" % d)
    PRT.pcol("[+] d2 = %d [mm]" % d2)
    PRT.pcol("[+] Mt = %d [mm]" % Mt)
    PRT.pcol("[+] l = %d [mm]" % l)

CLT = {
    "init": {
        "d": d,
        "d2": d2,
        "Mt": Mt,
        "l": l
    }
}

pi = math.pi

# Etapa [ I ]

# constante per etapa
micro_t = 0.08
niuA = 0.3
niuB = 0.3
EA = 2.1e5
EB = 2.1e5
E = 2.1e5
RmaxA = 5
RmaxB = 5
c = 1.1

# calcule folosind formule 
def calcKA(d,d1,d2,niuA):
    return (d**2 + d1**2)/(d**2 - d1**2) - niuA

def calcKB(d,d1,d2,niuB):
    return (d2**2 + d**2)/(d2**2 - d**2) + niuB

def calcPresiunea(Mt,micro_t, d, l):
    return (2*Mt*c) / (micro_t * pi * d**2 * l)

def calcSrelatiaLame(p_min, d, KA, KB, E):
    return p_min * d * 1e3 * (KA + KB)/E

def calcSminReal(s_min_nec, RmaxA, RmaxB):
    return s_min_nec + 2 * 0.6 * (RmaxA + RmaxB)

def verificareCondTransmSolic(X, AbatAlz, Arb, s_min_real, d):
    s_min_aj = X[AbatAlz]["date"][Arb]["ai"] - X[AbatAlz]["main"]["AS"]
    
    if s_min_aj > s_min_real:

        PRT.pcol("[!] Verificare pentru X = ø%d %s/%s" % (d,AbatAlz,Arb))
        PRT.pcol("[*] %s [AS = %d | AI = %d]" % (AbatAlz, X[AbatAlz]["main"]["AS"], X[AbatAlz]["main"]["AI"]))
        PRT.pcol("[*] %s [as = %d | ai = %d]" % (Arb, X[AbatAlz]["date"][Arb]["as"], X[AbatAlz]["date"][Arb]["ai"]))
        print("\n")
        PRT.pcol("[+] s min aj = %.3f [μm]" % s_min_aj)

        s_max_aj = X[AbatAlz]["date"][Arb]["as"] - X[AbatAlz]["main"]["AI"]
        PRT.pcol("[+] s max aj = %.3f [μm]" % s_max_aj)

        s_max_real = s_max_aj - 2 * 0.6 * (RmaxA + RmaxB)
        PRT.pcol("[+] s max real = %.3f [μm]" % s_max_real)

        return s_max_real
    else:
        # PRT.pcol("[-] Nu sa verificat conditia de transmitere a solicitarii Mt pentru X = ø%d %s/%s" % (d,AbatAlz,Arb))
        return None

def calcCoefSiguranta(s_max_real, d, d2, KA, KB):
    pmax = (s_max_real * 1e-3) / (d * (KA + KB)/E)
    PRT.pcol("[+] p max = %.3f [N/mm^2]" %  pmax)

    sigma_t = pmax * (d2**2 + d**2) / (d2**2 - d**2)
    PRT.pcol("[+] σt = %.3f [N/mm^2]" % sigma_t)

    sigma_r = -pmax
    PRT.pcol("[+] σr = %.3f [N/mm^2]" % sigma_r)

    sigma_ech = math.sqrt(sigma_t**2 + sigma_t * sigma_r + sigma_r**2)
    PRT.pcol("[+] σech = %.3f [N/mm^2]" % sigma_ech)
    return [sigma_ech, pmax]

def calcFortaAxiala(pmax, d, l):
    miu = 0.2 # 0.2 .. 0.25
    k = 1.1 # 1.1 .. 2

    N = pmax * pi * d * l
    Ffmax = miu * N
    FA = k * Ffmax

    PRT.pcol("[+] N = %.3f [N]" % N)
    PRT.pcol("[+] Ffmax = %.3f [N]" % Ffmax)
    PRT.pcol("[+] FA = %.3f [N]" % FA)
    FA = math.floor(FA)
    part = FA % 100
    if part <= 50:
        part = 50 - part
    else:
        part = 100 - part
    FAr = FA + part
    PRT.pcol("[+] FA rotunjit = %d [N]" % FAr)

    return FAr


# obtinere obiect date in baza de diametru
class AL:
    X = getDiametru(d)

PRT.header(":: ----- Etapa [I] ----- ::")
p = calcPresiunea(Mt, micro_t, d, l)
KA = calcKA(d, d1, d2, niuA)
KB = calcKB(d, d1, d2, niuB)
s = calcSrelatiaLame(p, d, KA, KB, E)

s_min_real = calcSminReal(s, RmaxA, RmaxB)

PRT.pcol("[+] KA = %.3f" % KA)
PRT.pcol("[+] KB = %.3f" % KB)
PRT.pcol("[+] p = %.3f" % p)
PRT.pcol("[+] s = %.3f [μm]" % s)
PRT.pcol("[+] s min real = %.3f [μm]" % s_min_real)

sigma_ech = {}

good = 0

AbatAlz = [key for key in AL.X]

for AbtAlz in AbatAlz:    
    for key in AL.X[AbtAlz]["date"]:
        if good < 4:
            s_max_real = verificareCondTransmSolic(AL.X, AbtAlz, key, s_min_real, d)
            if s_max_real is not None:
                sigma_ech[AbtAlz+"/"+key] = calcCoefSiguranta(s_max_real, d, d2, KA, KB)
                good+=1
        else: 
            break

coef_max = 0

tab_cf = []
corner = Fore.LIGHTYELLOW_EX+"Mat"+Style.RESET_ALL +"\\"+ Fore.LIGHTCYAN_EX+"Ajust"+Style.RESET_ALL
ls_ajustaje = [corner]

materiale = OLC.Marca

if good == 4:
    del materiale["OLC 15"]
    PRT.pcol("[!] Calculul coeficientului de siguranta pentru 3 marci de otel")
else:
    PRT.pcol("[!] Calculul coeficientului de siguranta pentru 4 marci de otel")

for marca in materiale:
    #PRT.pcol("[!] Pentru marca [ %s ] " % marca)
    
    line = [Fore.LIGHTYELLOW_EX+marca+Style.RESET_ALL]
    
    for ajustaj in sigma_ech:
        #PRT.pcol("[!] Pentru ajustajul [ %s ]  cu sigma_ech = %.3f" % (ajustaj, sigma_ech[ajustaj][0]))
        coef_sig = materiale[marca]["Rp02"]  / sigma_ech[ajustaj][0]
        #PRT.pcol("[*] Coef sig = %.3f" % coef_sig)
        
        if coef_max < coef_sig:
            coef_max = coef_sig
            material = marca
            ajs = ajustaj
            pmax = sigma_ech[ajustaj][1]
    
        if not ajustaj in ls_ajustaje:
            ls_ajustaje.append(Fore.LIGHTCYAN_EX+ajustaj+Style.RESET_ALL)

        line.append(coef_sig)

    tab_cf.append(line)
for tr in range(len(tab_cf)):
    for td in range(len(tab_cf[tr])):
        if coef_max == tab_cf[tr][td]:
            tab_cf[tr][td] = Fore.CYAN + str(tab_cf[tr][td]) + Style.RESET_ALL

print(tabulate(tab_cf,headers=ls_ajustaje))
PRT.pcol("[!] Concluzii: Materialul pentru arbore si butuc [%s], suprafata de asamblare este alezata, iar ajustajul este: ø%d %s" % (material, d, ajs))
FA = calcFortaAxiala(pmax, d, l)


"""----------------------------------------------------------------------------------------"""

# Etapa [ II ]

# constante
cf = 2.4 # 2.4 .. 4
c = 1.1  # 1.1 .. 2

miu_placuta = {"bronz": 0.08, "otel": 0.2}

# calcule folosind formule
def calc_d1(FA):
    return (56/pi) * math.sqrt((cf * FA) / (pi * E))

def calc_d1_stelat(FA, Rp02, fi):
    return math.sqrt((4*FA*c)/(pi*Rp02*fi))

def calc_alpha2(p,d2):
    return p/(pi * d2)

def calc_rho(miu, beta, filet):
    if filet == "patrat":
        return miu #math.atan(miu)
    else:
        return miu/0.965925826 #math.cos(beta/2)

def calc_psp(d0, t, di, FA):
    return FA/((pi/4) * ((d0 - 2*t)**2 - di**2)) 

def calc_sigma_asfr(Rp02, cstr):
    return Rp02/cstr

def calc_pmax(FA, R):
    return 0.39 * ((FA*E**2)/R**2)**(1/3)

def calc_r0(FA, R):
    return 1.11 * ((FA*R)/E)**(1/3)

# Momentele de torsiune
def calc_Mt2a(miu, psp, d0, t, di):
    return (pi/12) * miu * psp * ((d0 - 2*t)**3 - di**3)

def calc_Mt2b(miu, pmax, r0):
    return (pi**2 / 8) * miu * pmax * r0**3

def calc_Mt2c(FA, di):
    miu = 0.3 # 0.25 .. 0.3
    db = math.ceil(di * 1.5) # 1.25 .. 1.5 
    PRT.pcol("[+] Adoptam db = %.3f [mm]" % db)
    gama = 0.866025403
    return (miu * FA * db * gama)/4

def redimensionare(diam, filet):
    if filet == 'patrat':
        return getDiametru_d1_FiletPatrat(diam)
    else:
        return getDiametru_d1_FiletTrapez(diam)

def afisare_date_extrase_STAS(date, filet):
    PRT.pcol("[!] Date extrase din STAS pentru filetul cu profil [%s]:" % filet)
    PRT.pcol("[+] d1 = %.3f [mm]" % date[filet]["STAS"]["d1"])
    PRT.pcol("[+] d = %.3f [mm]" % date[filet]["STAS"]["d"])
    PRT.pcol("[+] d2 = D2 = %.3f [mm]" % date[filet]["STAS"]["d2"])
    PRT.pcol("[+] p = %.3f [-]" % date[filet]["STAS"]["p"])

def calculePerTipFilet(date, d1, FA, materiale):
    d1_max = d1
    # verificare pentru materiale diferite
    for key in materiale:
        PRT.pcol("[!] Verificare pentru [%s] :" % key)
        d1_stelat = calc_d1_stelat(FA, materiale[key]["Rp02"], materiale[key]["fi"])
        PRT.pcol("[+] d1* = %.3f [mm]" % d1_stelat)
        if d1_stelat > date["STAS"]["d1"]:
            PRT.pcol("[-] Conditia d1* > d1 STAS [Verificata] => Redemensionarea d1 STAS")
            if d1_max < d1_stelat:
                d1_max = d1_stelat
                material_selectat = key
        else:
            PRT.pcol("[+] Conditia d1* > d1 STAS [Nu este verificata] => ramane situatia initiala")
        
    return d1_max, material_selectat

def inc_placuta_miu(miu, placuta):
    
    if miu == 0:
        return miu_placuta[placuta]

    if placuta == 'bronz':
        if 0.08 < miu and miu+0.01 < 0.2:
            return miu + 0.01
        else:
            return miu 
    else:
        if 0.2 < miu and miu+0.01 < 0.25:
            return miu +0.01
        else:
            return miu 

def afisare_unghiuri(filet, alfa2, rho):
    PRT.pcol("[+] tg(α2) = %.8f" % alfa2)
    if filet == "patrat":
        PRT.pcol("[+] tg(ρ) = %.8f" % rho)
    else:
        PRT.pcol("[+] tg(ρ') = %.8f" % rho)

def verificareAutofranare(date, placuta, filet):
    PRT.pcol("[!] Verificare la autofranare pentru placuta de [%s] :" % placuta)
    placuta = placuta.split("/")[-1]
    alfa2 = calc_alpha2(date["STAS"]["p"], date["STAS"]["d2"])
    while True:
        miu = inc_placuta_miu(0,placuta)
        rho = calc_rho(miu, date["beta"], filet)

        if alfa2 <= rho:
            PRT.pcol("[+] Selectare: otel/%s miu = %.3f" % (placuta, miu))
            afisare_unghiuri(filet, alfa2, rho)
            PRT.pcol("[*] Verificare:  α2 <= ρ -> [%r]" % (alfa2<=rho))
            return miu

PRT.header(":: ----- Etapa [II] ----- ::")

# Determinarea tipodimensionala
d1 = calc_d1(FA)
PRT.pcol("[+] d1 calc = %.3f [mm]" % d1)

perFilet = {
    "patrat": {
        "STAS": getDiametru_d1_FiletPatrat(d1),
        "beta": 0,
    },

    "trapezoidal": {
        "STAS": getDiametru_d1_FiletTrapez(d1),
        "beta": pi/6,
    },
    
    "miu": {
        "otel/bronz": 0.08,
        "otel/otel": 0.2
    }
} 

afisare_date_extrase_STAS(perFilet, "patrat")
afisare_date_extrase_STAS(perFilet, "trapezoidal")

materiale = {
    "OL37":  {"fi": 0.75, "Rp02": 210},
    "OLC45": {"fi": 0.89, "Rp02": 360}
}

PRT.header(":: Determinarea tipodimensionala ::")
d1_verificare, material_selectat = calculePerTipFilet(perFilet["patrat"], d1, FA, materiale)
if d1_verificare != perFilet["patrat"]["STAS"]["d1"]:
    for filet in ["patrat", "trapezoidal"]:
        perFilet[filet]["STAS"] = redimensionare(d1_verificare, filet)

        PRT.pcol("[!] Redemensionare d1 STAS filet %s " % filet)
        afisare_date_extrase_STAS(perFilet, filet)


for key in ["patrat", "trapezoidal"]:
    PRT.header(":: Calcule pentru profil [%s] ::" % key)
    # verificare la autofranare
    for placuta in perFilet["miu"]:
        perFilet["miu"][placuta] = verificareAutofranare(perFilet[key], placuta, key)

# calculul momentului de frecare din pivot
PRT.pcol("[!] Calculul momentului de frecare din pivot")


csfr = 1.5 # 1.5 .. 1.8
Rp02_bz = 300
d0 = perFilet["patrat"]["STAS"]["d1"]
di = getDiametru_di(d0)
t = 1 # 1 .. 2
psp = math.ceil(calc_psp(d0, t, di, FA))

sigma_asfr = math.trunc(calc_sigma_asfr(Rp02_bz, csfr))

PRT.pcol("[+] csfr = %.3f" % csfr)
PRT.pcol("[+] Rp02[bz] = %.3f [MPa] -> bronz aliat" % Rp02_bz)
PRT.pcol("[+] d0 = [d1] = %.3f [mm]" % d0)
PRT.pcol("[+] di = %.3f [mm]" % di)
PRT.pcol("[+] t = %.3f" % t)
PRT.pcol("[+] psp = %.3f" % psp)
PRT.pcol("[+] σsfr = %.3f" % sigma_asfr)
PRT.pcol("[*] Verificare psp <= σsfr : [%r] " % (psp <= sigma_asfr))

Mt2a = calc_Mt2a(perFilet["miu"]["otel/bronz"], psp, d0, t, di)
PRT.pcol("[+] Mt2a = %.3f " % Mt2a)

# sprijin pe placute din otel cu suprafata sferica
PRT.pcol("[!] Sprijin pe placute din otel cu suprafata sferica")
R = 100 # 100 .. 200 
pmax = calc_pmax(FA, R)
r0 = calc_r0(FA, R)
PRT.pcol("[+] R = %.3f" % R)
PRT.pcol("[+] pmax = %.3f" % pmax)
PRT.pcol("[+] r0 = %.3f" % r0)

Mt2b = calc_Mt2b(perFilet["miu"]["otel/otel"], pmax, r0)
PRT.pcol("[+] Mt2b = %.3f " % Mt2b)

# sprijin prin intermediul unei bile din otel 
PRT.pcol("[!] Sprijin prin intermediul unei bile din otel")

Mt2c = calc_Mt2c(FA, di)
PRT.pcol("[+] Mt2c = %.3f " % Mt2c)

perFilet["momente"] = {
    "Mt2a": Mt2a,
    "Mt2b": Mt2b,
    "Mt2c": Mt2c
}

"""-------------------------------------------------------------------------------"""

# Etapa [ III ]

def suma_unghiuri(alfa2, rho):
    return (alfa2+rho) / (1- (alfa2*rho))

def calc_tau_torsiune(Mt12,d1):
    return (16 * Mt12) / (pi * d1**3)

def calc_tau_adm_tors(Rp02):
    return Rp02/2

def calc_sigma_c(FA, d1_STAS):
    return (4*FA) / (pi * d1_STAS**2)

def calc_sigma_ech(sigma_c, tau_t):
    return math.sqrt( sigma_c**2 + 3 * tau_t**2 )


# Etapa [ III ]

PRT.header(":: ----- Etapa [III] ----- ::")

profile = ["patrat", "trapezoidal"]
placute = ["bronz", "otel"]

aps = ["'", "''"]
ind = ["a","b"]

PRT.header(":: Alegerea variantei optime de valoarea randamentului ::")

dateFilete = perFilet

dateFilete["patrat"]["momente1"] = {}
dateFilete["trapezoidal"]["momente1"] = {}

randamente = []
max_rand = 0

for i in range(2):
    PRT.header(":: Pentru contact [otel-%s] ::" % placute[i])

    for j in range(2):
        profil_actual = profile[j]
        PRT.pcol("[!] La profil [%s] : d2 = %.2f [mm]  FA = %d [N] p = %d [mm] ->" % (profil_actual, dateFilete[profil_actual]["STAS"]["d2"], FA, dateFilete[profil_actual]["STAS"]["p"]))
        tag_alpha = calc_alpha2(dateFilete[profil_actual]["STAS"]["p"], dateFilete[profil_actual]["STAS"]["d2"])
        tag_rho = calc_rho(dateFilete["miu"]["otel/%s" % placute[i]], dateFilete[profil_actual]["beta"], profil_actual)
        
        sm_ar = suma_unghiuri(tag_alpha, tag_rho)

        Mt_xy = FA * (dateFilete[profil_actual]["STAS"]["d2"]/2) * sm_ar

        dateFilete[profil_actual]["momente1"]["Mt%s1%s" % (aps[i], ind[j])] = Mt_xy
        
        afisare_unghiuri(profil_actual, tag_alpha, tag_rho)
        PRT.pcol("[+] tg(α2 + ρ) = %.8f" % sm_ar)
        PRT.pcol("[+] Mt%s1%s = %.3f [N*mm]" % (aps[i], ind[j], Mt_xy))

        for Mt2 in dateFilete["momente"]:
            Mt2_temp = dateFilete["momente"][Mt2]
            rand = 100 * (FA * dateFilete[profil_actual]["STAS"]["d2"] * tag_alpha) / (2*(Mt_xy + Mt2_temp))
            randamente.append(rand)

            if rand > max_rand:
                max_rand = rand
                Mt1x = Mt_xy
                Mt2y = Mt2_temp
                tp_filet = profil_actual
                momente = Mt2                    
                mater_placuta = placute[i]

PRT.pcol("[!] Calcul randamente :")
for i in range(len(randamente)):
    PRT.pcol("[+] η%d = %.2f %%" % (i+1, randamente[i]))

d1_STAS = dateFilete[tp_filet]["STAS"]["d1"]

material_surub = material

PRT.header(":: Concluzii ::")
PRT.pcol("- materialul surubului : otel -> "+material_surub)
PRT.pcol("- materialul piulitei : "+mater_placuta)
PRT.pcol("- tipul filetului si dimensiunea : filet %s cu d1 = %.2f [mm]" % (tp_filet, d1_STAS))
PRT.pcol("- tipul sprijinului : " + sprijin(momente))

Rp02 = OLC.Marca[material_surub]["Rp02"]

PRT.header(":: Calculul de rezistenta a surubului ::")
PRT.pcol("a) Rp02 = %d [N/mm^2]; d1 (STAS) = %.2f [mm]" % (Rp02, d1_STAS))
PRT.pcol("[+] η max = %.2f %%" % max_rand)
PRT.pcol(f"[+] Mt1x = {Mt1x} | Mt2y = {Mt2y}")
tau_t = calc_tau_torsiune(Mt1x+Mt2y, d1_STAS)
tau_at = calc_tau_adm_tors(Rp02)
PRT.pcol("[+] τt = %.3f [MPa]" % tau_t)
PRT.pcol("[+] τa,t = %.3f [MPa]" % tau_at)
if tau_t <= tau_at:
    PRT.pcol("[*] τt <= τa,t -> [verifica]")

PRT.pcol("b) FA = %d [N]; Rp02 = %d [N/mm^2]; d1 (STAS) = %.2f [mm]; τa,t = %.3f [MPa]" % (FA, Rp02, d1_STAS, tau_at))
sigma_c = calc_sigma_c(FA, d1_STAS)
tau_t_b = calc_tau_torsiune(Mt2y, d1_STAS)
sigma_ech = calc_sigma_ech(sigma_c, tau_t_b)


PRT.pcol("[+] σc = %.3f [MPa]" % sigma_c)
PRT.pcol("[+] τt = %.3f [MPa]" % tau_t_b)
PRT.pcol("[+] σech = %.3f [MPa]" % sigma_ech)
if sigma_ech <= tau_at:
    PRT.pcol("[*] σech <= σa -> [verifica]")
else:
    PRT.pcol("[-] σech <= σa -> [nu verifica]")

"""-----------------------------------------------------------------------------------"""

# Etapa [ IV ]

def verificari_tensiuni(tens_ef, den_tens_ef, tens_a, den_tens_a):
    PRT.pcol("[+] %s = %.3f [Mpa]" % (den_tens_ef, tens_ef))
    PRT.pcol("[+] %s = %.3f [Mpa]" % (den_tens_a, tens_a))

    if tens_ef <= tens_a:
        PRT.pcol("[*] Sa verificat conditia %s <= %s" % (den_tens_ef, den_tens_a))
        return True
    else:
        PRT.pcol("[-] Nu sa verificat conditia %s <= %s" % (den_tens_ef, den_tens_a))
        return False 

PRT.header(":: ----- Etapa [IV] ----- ::")

datFil = dateFilete["patrat"]["STAS"]
datFil["D1"] = datFil["d1"]

PRT.header(":: 1. Calculul numarului de spire ::")

D = datFil["d"] + 0.5

PRT.pcol("[+] D = %.1f [mm]; D1 = d1 = %.2f [mm]; p = %d [mm]; " % (D,datFil["d1"],datFil["p"]))
PRT.pcol("[+] d = %d [mm]; FA = %d [N]; " % (datFil["d"],FA))

a = 1
csi = 1.8 # 1.5 .. 1.8

Rp02 = Rp02_bz
z = (4*csi*FA) / (Rp02 * pi * (datFil["d"]**2 - datFil["D1"]**2))

PRT.pcol("[!] Pentru piulita se alege [ %s ] cu Rp02 = %d " % (mater_placuta,Rp02))
PRT.pcol("[*] z = %.3f => z min = %d => admitem z = %d" % (z,math.ceil(z),6))
z_min = math.ceil(z)
z = 8

De = math.sqrt((8 * FA + Rp02 * pi * datFil["d1"]**2) / (pi * Rp02))

PRT.pcol("[+] De = %.2f [mm] " % De)
if De < d + 12:
    De = d+12 +1
De = math.ceil(De)

PRT.pcol("[*] Adoptam De = %d [mm]" % De)

D0 = 2 * De # *2 .. 2.5
h = z * datFil["p"]; m = 0.25 * h; r = 2; t = 1

PRT.pcol("[+] D0 = %.2f [mm]; m = %.2f [mm]; h = %.2f [mm]; r = %d [mm]; t = %d [mm]" % (D0, m,h,r,t))

n = 2

Mt_stelat = Mt1x + Mt2y
PRT.pcol(f"[*] Mt* = {Mt_stelat} [N*mm]")

Ds = De + (D0-De)/2
ds = math.sqrt((16 * Mt_stelat)/(0.65 * Ds * Rp02 * n * pi))
PRT.pcol("[+] Ds = %.2f [mm]" % Ds)
PRT.pcol("[+] ds = %.2f [mm]" % ds)

ds = math.ceil(ds)
PRT.pcol("[*] ds adoptat = %d [mm]" % ds)
PRT.pcol("[+] n = %d [ nr suruburi]" % n)

if datFil["d"] + 12 > De:
    De_desen = datFil["d"] + 12 + z_min
else:
    De_desen = De
PRT.pcol("[!] Se adopta pentru desen: De = %d [mm]; h = %d [mm];" % (De_desen, h))

PRT.header(":: Verificari din conditia de rezistenta pentru filetul si corpul piulitei ::")
PRT.pcol("a) Verificarea spirei filetului la incovoiere :")
sigma_ef_inc = (24*FA*((datFil["d"]-datFil["D1"])/4 + a)) / (pi * D * h * datFil["p"])
csi = 1.5
PRT.pcol("[*] c s,i = %.1f" % csi)
sigma_a_inc = Rp02/csi

if not verificari_tensiuni(sigma_ef_inc,"σef inc",sigma_a_inc,"σa inc"):
    while True:
        h += 1
        sigma_ef_inc = (24*FA*((datFil["d"]-datFil["D1"])/4 + a)) / (pi * D * h * datFil["p"])
        if sigma_ef_inc <= sigma_a_inc:
            PRT.pcol(f"[!] Redimensionare h = {h}")
            break

PRT.pcol("b) Verificarea spirei filetului la forfecare :")
tau_ef_f = (2/pi)*(FA/(datFil["d"]*h))
tau_a_f = 0.65 * Rp02 / 1.5

verificari_tensiuni(tau_ef_f,"τef f",tau_a_f,"τa f")

PRT.pcol("c) Verificarea gulerului de reazem la strivire :")
sigma_ef_str = FA / ((pi/4) * ((D0-2*t)**2 - (De_desen+2*r)**2)) 
sigma_a_str = Rp02/1.5

verificari_tensiuni(sigma_ef_str,"σef str",sigma_a_str,"σa str")

PRT.pcol("d) Verificarea gulerului de reazem la efort compus :")
sigma_ef_inc = (3 * FA * (D0-De)) / (2 * pi * De * m**2)
tau_ef_f = FA / (pi * De * m)
PRT.pcol(f"[+] Sigma ef inc = {sigma_ef_inc} [MPa];  tau ef f = {tau_ef_f} [MPa]")
sigma_ech_ef = math.sqrt(sigma_ef_inc**2 + 3 * tau_ef_f**2)
sigma_a_t = Rp02/2

if not verificari_tensiuni(sigma_ech_ef, "σech ef", sigma_a_t, "σa t"):

    while not sigma_ech_ef < sigma_a_t:
        m += 1
        sigma_ef_inc = (3 * FA * (D0-De)) / (2 * pi * De * m**2)
        tau_ef_f = FA / (pi * De * m)
        sigma_ech_ef = math.sqrt(sigma_ef_inc**2 + 3 * tau_ef_f**2)
        sigma_a_t = Rp02/2

    
    PRT.pcol(f"[!] >>> Redimensionare m = {m}")
    PRT.pcol(f"[+] Sigma ef inc = {sigma_ef_inc} [MPa];  tau ef f = {tau_ef_f} [MPa]")
    verificari_tensiuni(sigma_ech_ef, "σech ef", sigma_a_t, "σa t")


# Etapa [ V ]
PRT.header(":: ----- Etapa [V] ----- ::")

PRT.header(":: Proiectarea boltului ::")

PRT.pcol("- De (diametrul corpului piulitei) = %d [mm];" % De)
PRT.pcol("- D0 (diametrul gulerului piulitei) = %d [mm];" % D0)
PRT.pcol("- db (diametrul boltului)")

PRT.pcol("[!] Se adopta pentru bolt, brat extractor si cadru un [OLC 45]")
PRT.pcol("a) Calcularea diametrului boltului")

Rp02 = 360
PRT.pcol("[+] FA = %d [N]; D0 = %d [mm]; Rp0,2 = %d [N/mm^2]" % (FA, D0, Rp02))

db1 = ((1.8 * 32 * D0 * FA)/(pi * Rp02 * 12))**(1/3)
db2 = math.sqrt((2*FA)/(0.65*Rp02*pi))
db3 = (2*FA)/(pi*D0*Rp02)

PRT.pcol("[+] db.1 = %.2f [mm]" % db1)
PRT.pcol("[+] db.2 = %.2f [mm]" % db2)
PRT.pcol("[+] db.3 = %.2f [mm]" % db3)

db = max(db1,db2,db3)
PRT.pcol("[+] db = max(db.1; db.2; db.3) => db = %.2f [mm]" % db)
if db < math.ceil(db):
    db = math.ceil(db)  
    PRT.pcol("[!] Se majoreaza in adaos valoarea diametrului boltului, deci db = %d [mm]" % db)

PRT.pcol("b) Schita boltului cu dimensiunile calculate si adoptate")

PRT.pcol("[*] t = (1..2) [mm]; ds = (2..4) [mm]; r = (2..3) [mm]; D0 = %d [mm]; db = %d [mm];" % (D0,db))

ltb = D0 + 10
lcb = db
dcb = 2 * db

PRT.pcol("[+] ltb = %d [mm]" % (ltb))
PRT.pcol("[+] lcb = db = %d [mm]; dcb = %d [mm];" % (lcb,dcb))

PRT.header(":: Proiectarea bratelor extractoare ::")
PRT.pcol("[*] D0 = %d [mm]; Rp0,2 = %d [N/mm^2]; FA = %d [N]" % (D0, Rp02,FA))

x = (3*FA) / (2*D0*Rp02)
PRT.pcol("[+] x = %.2f [mm]" % x)

if x < math.ceil(x):
    x = math.ceil(x)
    msg = " si se va majora la o valoare intreaga, adica x = %d [mm]" % x
else:
    msg = ""

PRT.pcol("[!] Aceasta valoarea x este valoarea minima necesara pentru ca sa nu se rupa bratul extractor in dreptul gaurii" + msg)

y1 = 2 # 2..5
y2 = 2 # 2..5
PRT.pcol("[*] y1 = %.2f & y2 = %.2f [mm]; d2 = %d [mm]; d = %d [mm]" % (y1,y2, CLT["init"]["d2"], d))

g1 = ((9 * 1.8 * FA * (CLT["init"]["d2"] - CLT["init"]["d"] - 2*y2 + 2*y1))/(2 * D0 * Rp02))**(1/2)
g2 = (3*FA)/(0.65 * D0 * Rp02)
g = max(g1,g2)

PRT.pcol("[+] g1 = %.2f [mm]" % g1)
PRT.pcol("[+] g2 = %.2f [mm]" % g2)
PRT.pcol("[+] g = %.2f [mm] => g = %d [mm]" % (g,math.ceil(g)))
g = math.ceil(g)

PRT.pcol("[!] Se face verificarea")
sigma_str = (3*FA) /(D0* (CLT["init"]["d2"] - CLT["init"]["d"] - 2*y2))
sigma_a_str = Rp02/1.8

verificari_tensiuni(sigma_str,"σ str",sigma_a_str,"σ a str")

PRT.header(":: ----- Etapa [VI] ----- ::")

PRT.header(":: Proiectarea cadrului ::")
PRT.pcol("[!] Se adopta urmatoarele dimensiuni in [mm] : ")

a = db/2 + x + y1 + CLT["init"]["d2"]/2

while True:
    sigma_a_i = Rp02/1.8
    sigma_i1 = (3 * FA * a)/((D0-De)*h**2)
    if sigma_i1 <= sigma_a_i:
        break
    else:
        h += 1

PRT.pcol("[*] D0 = %.2f; De = %.2f; m = %.2f; g = %.2f; db = %.2f;\n    h = %.2f; x = %d; a = %.2f; d2 = %.2f." % (D0, De, m, g, db, h, x, a, CLT["init"]["d2"]))
print("\n")

if "1" == INP.icol("Pozitionare cadru 1:(h) / 2:(h-m) -> "):
    hp = h
    varianta_cadru = "(h)"
else:
    hp = h-m
    varianta_cadru = "(h-m)"

PRT.pcol("[!] (1..1) solicitare de incovoiere")
sigma_a_i = Rp02/1.8
sigma_i1 = (3 * FA * a)/((D0-De)*h**2)
verificari_tensiuni(sigma_i1,"σ i1", sigma_a_i,"σ a i")

PRT.pcol("[!] (2..2) solicitare de incovoiere")
sigma_i2 = (3 * FA * (2*a-De))/(D0*(h-m)**2)
verificari_tensiuni(sigma_i2,"σ i2", sigma_a_i,"σ a i")

PRT.pcol("[!] (3..3) solicitare compusa de incovoiere si forfecare")
sigma_i3 = absolute((3 * FA * (2*a-D0))/(D0*(h-m)**2))
tau_f3 = FA/(2*D0*(h-m))
PRT.pcol("[+] σ i3 = %.3f [MPa]" % sigma_i3)
PRT.pcol("[+] τ f3 = %.3f [MPa]" % tau_f3)

sigma_ech = math.sqrt(sigma_i3**2 + 3*tau_f3**2)
sigma_a_t = Rp02/2
verificari_tensiuni(sigma_ech,"σ ech",sigma_a_t, "σ a t")

PRT.pcol("[!] (4..4) solicitare de incovoiere")
sigma_i4 = (3*FA*(db/2 + x + y1))/((D0-D0/3) * h**2)
verificari_tensiuni(sigma_i4,"σ i4",sigma_a_i,"σ a i")

sigma_str = (3*FA)/(2*2*pi*db*D0)
verificari_tensiuni(sigma_str,"σ str",Rp02/1.8,"Rp0,2/1.8")

PRT.header(":: ----- Etapa [VII] ----- ::")

PRT.header(":: Proiectarea formei constructive a surubului ::")

PRT.pcol("- lca - lungime cap de actionare")
PRT.pcol("- li - lungime tronson intermediar ")
PRT.pcol("- lf - lungime filetata")
PRT.pcol("- ls - lungime tronson sprijin")


PRT.pcol(f"[*] Mt1x + Mt2y => {Mt1x} + {Mt2y} = {Mt1x+Mt2y}")
MT_xy = math.ceil((Mt1x + Mt2y)/1000)*1000
PRT.pcol(f"[*] Mt = {MT_xy}")
b = 1000 # [mm] max value 1000
F = math.ceil(MT_xy) / b

PRT.pcol("[+] b = %d [mm]" % b)
PRT.pcol("[+] F = %.3f [mm]" % F)

d_stelat = datFil["d1"]
PRT.pcol("[+] d* = [d1 filet] = %.2f [mm]" % datFil["d1"])

PRT.pcol("[!] Situatii de actionare :")

Mtc = Mt_xy

# hexagon
l6 = d_stelat/2 
b6 = l6/3

PRT.pcol("[+] l6 = %.2f; b6 = %.2f [mm]" % (l6,b6))

# patrat
l4 = d_stelat/math.sqrt(2)
b4 = l4/3
PRT.pcol("[+] l4 = %.2f; b4 = %.2f [mm]" % (l4,b4))

# triunghi
l3 = (d_stelat*3)/(2*math.sqrt(3))
b3 = l3/3
PRT.pcol("[+] l3 = %.2f; b3 = %.2f [mm]" % (l3,b3))

def calc_afisare_F_h(cheia,sectiunea,x,Mtc,bx,lx,Rp02):
    PRT.pcol("[!] Cheie %s : sectiunea %s [h%s]" % (cheia,sectiunea,x))
    Fx = Mtc/bx
    hx = (2*Fx)/(lx*Rp02)
    PRT.pcol("[+] F%s = %.2f [N]" % (x,Fx))
    PRT.pcol("[+] h%s = %.2f [mm]" % (x,hx))
    return Fx, hx

F6, h6 = calc_afisare_F_h("fixa","hexagonala","6",Mtc,b6,l6,Rp02)
F4, h4 = calc_afisare_F_h("fixa","patrata","4",Mtc,b4,l4,Rp02)

F6_stelat, h6_stelat = calc_afisare_F_h("inelara","hexagonala","6*",Mtc,3*b6,l6,Rp02)
F4_stelat, h4_stelat = calc_afisare_F_h("inelara","patrata","4*",Mtc,3*b4,l4,Rp02)
F6_stelat, h6_stelat = calc_afisare_F_h("inelara","triunghiulara","3*",Mtc,3*b3,l3,Rp02)

Hmin = g + CLT["init"]["l"] + 20 + 15 + hp

PRT.header(":: Concluzii ::")
PRT.pcol("[!] Surub: %s, d1 = %.2f mm, d2 = %.2f mm, d = %d mm, p = %d mm" % (tp_filet,datFil["d1"], datFil["d2"], datFil["d"], datFil["p"]))
PRT.pcol("[!] Piulita: %s, D2=d2= %.2f mm, z min = %d mm, p = %d mm, D0 = %d mm, De = %d mm, m = %d mm" % (tp_filet, datFil["d2"],z,datFil["p"],D0,De_desen,m))
PRT.pcol("[!] Bolt: db = %d mm, dcb = %.2f mm, ltb = %.2f mm" % (db,dcb,ltb))
PRT.pcol("[!] Brat extractor: g = %d mm, x = %d mm, Hmin = %d mm" % (g,x,Hmin))
PRT.pcol("[!] Cadru: D0 = %d mm, a = %.2f mm, h_cadru = %d mm in varianta %s" % (D0,a,h,varianta_cadru))