# coding: utf-8
import numpy as np
from lxml import etree
from classes import OpticalElement, Beamline

tree = etree.parse(r"D:\Dennetiere\optix\bin\test\system.xml")
beamline = Beamline()
for user in tree.xpath("/system/element"):
    new_element = OpticalElement(name=user.get("name"), next=user.get("next"), previous=user.get("previous"))
    beamline.add_element(new_element)
    for param in user.xpath("/"):
        print(param.text)
for element in beamline.elements:
    print(element)
print("#"*80)
beamline.chain()

for element in beamline.elements:
    print(element)
print("#"*80)
for chain in beamline.chains:
    desc = ""
    for element in chain:
        desc += element.name+" -> "
    print(desc)