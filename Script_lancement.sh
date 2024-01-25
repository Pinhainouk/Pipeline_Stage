#!/bin/bash

fichier_liste="/home/elodie/Documents/Elodie/Config.txt" # Fichier contenant la liste des echantillons
script="/home/elodie/Documents/Elodie/Pipeline_Stage/Script.sh"

while IFS= read -r ligne
do
  # Vérifier si la ligne n'est pas vide
    if [ -n "$ligne" ]; then
        bash "$script" -o "/home/elodie/Documents/" -s "$ligne"
    fi
done < "$fichier_liste"
