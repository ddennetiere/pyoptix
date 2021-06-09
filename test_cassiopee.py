# coding: utf-8
import numpy as np
from lxml import etree

tree = etree.parse(r"D:\Dennetiere\optix\bin\test\system.xml")
for user in tree.xpath("/system/element"):
    print(user.get("name"))
    print(user.get("class"))
    print(user.get("next"))
    for param in user.xpath("/"):
        print(param.text)