#!/bin/bash
#
#   Script that parses EGAZ-XML files to extract the declared genome
#
#   Last Modified; June/30/02021
#
#   Version 0.0.1
#
#   Copyright (C) 2021 Manuel Rueda (manuel.rueda@crg.eu)
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, see <https://www.gnu.org/licenses/>.
#
#   If this program helps you in your research, please cite

set -eu

for xml in EGAZ/*xml
do 
 base=$(basename $xml .xml)
 echo -n "$base "
 perl parse_egaz_xml.pl $xml
done

# ================
#     Results
# ================
#   88687 (77.9 %) GRCh37
#   16075 (14.1 %) NA
#    8421 ( 7.4 %) hs37d5
#     669 ( 0.6 %) GRCh38
# ================
# Total 113852
