# -*- coding: utf-8 -*-

import csv, math, copy
from openbabel import *

conv = OBConversion()
conv.SetInAndOutFormats("smi", "mol2")

builder = OBBuilder()

ff = OBForceField.FindForceField("uff") # Fails with mmff94s :(
ff.SetLogLevel(OBFF_LOGLVL_NONE) # OBFF_LOGLVL_LOW
ff.SetLogToStdErr()

elements = OBElementTable()

def make3d(smiles):
    # SMILES -> 3D + optimization + charges + centering
    mol  = OBMol()
    conv.ReadString(mol, smiles)
    mol.AddHydrogens()
    builder.Build(mol)
    ff.Setup(mol)
    # ff.ConjugateGradients(200, 0.001)
    ff.SteepestDescent(50, 0.01)
    # Добавить конформационный анализ?
    # Сохранять в базе?
    ff.GetCoordinates(mol)
    ff.GetPartialCharges(mol)
    mol.Center()
    mol.SetEnergy(ff.Energy())
    return mol

def distance(org_atom, tag_atom):
    # Для удобства: расстояние между атомами
    ps = org_atom.x(), org_atom.y(), org_atom.z()
    qs = tag_atom.x(), tag_atom.y(), tag_atom.z()
    s  = sum([(p - q)**2 for p, q in zip(ps, qs)])
    return math.sqrt(s)

def dgeom(mol):
    # Дистанционная геометрия: расстояния
    atoms = [a for a in OBMolAtomIter(mol)]
    ds    = []
    for i in range(len(atoms) - 1):
        for j in range(i + 1, len(atoms)):
            ds.append(distance(atoms[i], atoms[j]))
    return ds

def box(mol):
    # Объем и всё, что около

    # В какой бокс можно вписать молекулу?
    c  = 1.0 # "Запас", Å
    xs = [a.x() for a in OBMolAtomIter(mol)]
    ys = [a.y() for a in OBMolAtomIter(mol)]
    zs = [a.z() for a in OBMolAtomIter(mol)]
    max_x = max(xs) + c
    max_y = max(ys) + c
    max_z = max(zs) + c
    min_x = min(xs) - c
    min_y = min(ys) - c
    min_z = min(zs) - c
    lx = max_x - min_x
    ly = max_y - min_y
    lz = max_z - min_z
    box_volume = lx * ly * lz

    # Заполнение бокса атомами, на сетке
    def t(x, y, z): # Тест на принадлежность узла атому из молекулы
        for a in OBMolAtomIter(mol):
            d_square = (x-a.x())**2 + (y-a.y())**2 + (z-a.z())**2
            r_square = elements.GetVdwRad(a.GetAtomicNum()) ** 2
            if d_square <= r_square:
                return True
        return False

    dl = 0.25   # Шаг сетки, Å
    nt = 0      # Узлов всего (nodes, total)
    nh = 0      # Узлов, принадлежащих атомам (nodes, hidden)

    x = min_x
    while x <= max_x:
        y = min_y
        while y <= max_y:
            z = min_z
            while z <= max_z:
                if t(x, y, z): nh += 1
                nt += 1
                z  += dl
            y += dl
        x += dl

    # То, что можно получить прямо из OBMol
    mol_wt = mol.GetMolWt()

    return (
        mol_wt,
        box_volume,
        mol_wt/box_volume,
        box_volume * nh/nt,
        float(nh)/nt)

def from_center(mol):
    # Для центрированной молекулы: расстояния от атомов до центра
    return [
         math.sqrt(a.x()**2 + a.y()**2 + a.z()**2) for a in OBMolAtomIter(mol)]

def center_of_mass(mol):
    # Для центрированной молекулы: расстояние центра масс от O

    def to_center_of(tuples):
        # Центр масс
        if not tuples:
            return 0.0, 0.0, 0.0
        ws = [w for _, _, _, w in tuples]
        sw = sum(ws)
        cx = sum([x * w for x, _, _, w in tuples]) / sw
        cy = sum([y * w for _, y, _, w in tuples]) / sw
        cz = sum([z * w for _, _, z, w in tuples]) / sw
        return cx, cy, cz # math.sqrt(cx**2 + cy**2 + cz**2)

    weighted = [(a.x(), a.y(), a.z(), a.GetAtomicMass())
        for a in OBMolAtomIter(mol)]

    min_m = min([m for _, _, _, m in weighted]) # Usually H
    max_m = max([m for _, _, _, m in weighted])

    cmx, cmy, cmz = to_center_of(weighted)
    dc_m = math.sqrt(cmx **2 + cmy **2 + cmz **2)

    return min_m, max_m, dc_m

def dipole(mol):
    # Молекула как диполь
    charged = [(a.x(), a.y(), a.z(), a.GetPartialCharge())
        for a in OBMolAtomIter(mol)]

    dip_x = sum([x * q for x, _, _, q in charged])
    dip_y = sum([y * q for _, y, _, q in charged])
    dip_z = sum([z * q for _, _, z, q in charged])
    dip_l = math.sqrt(dip_x**2 + dip_y**2 + dip_z**2)

    charged.sort(lambda a, b: cmp(a[3], b[3]))

    # Вектор "максимальный положительный заряд"
    pos_x, pos_y, pos_z, q = charged[len(charged)-1]
    dip_pos_l  = math.sqrt(pos_x**2 + pos_y**2 + pos_z**2)
    dip_pos_lq = math.sqrt(
        (q * pos_x)**2 +
        (q * pos_y)**2 +
        (q * pos_z)**2 )

    # Вектор "максимальный (по модулю) отрицательный заряд"
    neg_x, neg_y, neg_z = 0.0, 0.0, 0.0
    dip_neg_l  = 0.0
    dip_neg_lq = 0.0
    if charged[0][3] < 0.0:
        neg_x, neg_y, neg_z, q = charged[0]
        dip_neg_l  = math.sqrt(neg_x**2 + neg_y**2 + neg_z**2)
        dip_neg_lq = math.sqrt(
            (abs(q) * neg_x)**2 +
            (abs(q) * neg_y)**2 +
            (abs(q) * neg_z)**2 )

    # Угол между векторами
    # "максимальный положительный заряд" --
    # "максимальный (по модулю) отрицательный заряд"
    dip_angle = 0.0
    if dip_pos_l and dip_neg_l:
        prod = pos_x*neg_x + pos_y*neg_y + pos_z*neg_z
        try:
            dip_angle = math.acos(prod/(dip_pos_l * dip_neg_l))
        except ValueError:
            return 0.0

    return dip_l, dip_pos_l, dip_pos_lq, dip_neg_l, dip_neg_lq, dip_angle

def surface(mol):
    # Поверхность на сетке (медленно, неточно...)
    pass

# Вспомогательные статистические и прочие функции

def median(xs):
    ys = sorted(xs)
    n  = len(ys)
    i  = n // 2
    if n % 2: # Нечетное число элементов
        return ys[i]
    else: # Четное число элементов
        return (ys[i-1] + ys[i])/2

def average(xs):
    return sum(xs)/len(xs)

def rms(xs):
    return math.sqrt(sum([x**2 for x in xs])/len(xs))

# Hасчет дескрипторов

def count_three_dimensional(smiles):

    arr_3d_descriptors = {}

    # SMILES -> 3D + optimization
    mol = make3d(smiles)

    # Descriptors, distance geometry
    ds    = dgeom(mol)
    arr_3d_descriptors['max_d'] = max(ds)
    arr_3d_descriptors['min_d'] = min(ds)
    arr_3d_descriptors['avg_d'] = average(ds)
    arr_3d_descriptors['med_d'] = median(ds)
    arr_3d_descriptors['rms_d'] = rms(ds)

    # Descriptors, distances from center
    vs    = from_center(mol)
    arr_3d_descriptors['max_v'] = max(vs)
    arr_3d_descriptors['min_v'] = min(vs)
    arr_3d_descriptors['avg_v'] = average(vs)
    arr_3d_descriptors['med_v'] = median(vs)
    arr_3d_descriptors['rms_v'] = rms(vs)

    # Descriptors, partial charges
    ps    = [a.GetPartialCharge() for a in OBMolAtomIter(mol)]
    arr_3d_descriptors['max_p'] = max(ps)
    arr_3d_descriptors['min_p'] = min(ps)

    # Descriptors, volume
    mol_wt, box_volume, wt_vol_ratio, mol_volume, volume_ratio = box(mol)
    arr_3d_descriptors['mol_wt'] = mol_wt
    arr_3d_descriptors['box_volume'] = box_volume
    arr_3d_descriptors['wt_vol_ratio'] = wt_vol_ratio
    arr_3d_descriptors['mol_volume'] = mol_volume
    arr_3d_descriptors['volume_ratio'] = volume_ratio

    # Descriptors, center of mass
    min_m, max_m, dc_m = center_of_mass(mol)
    arr_3d_descriptors['min_m'] = min_m
    arr_3d_descriptors['max_m'] = max_m
    arr_3d_descriptors['dc_m'] = dc_m

    # Descriptors, dipole
    dip_l, dip_pos_l, dip_pos_lq, dip_neg_l, dip_neg_lq, dip_angle = dipole(mol)
    arr_3d_descriptors['dip_l'] = dip_l
    arr_3d_descriptors['dip_pos_l'] = dip_pos_l
    arr_3d_descriptors['dip_pos_lq'] = dip_pos_lq
    arr_3d_descriptors['dip_neg_l'] = dip_neg_l
    arr_3d_descriptors['dip_neg_lq'] = dip_neg_lq
    arr_3d_descriptors['dip_angle'] = dip_angle

    # Return array of descriptors
    return arr_3d_descriptors
