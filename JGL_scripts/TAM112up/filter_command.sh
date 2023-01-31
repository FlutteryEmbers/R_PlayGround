

### awk -F "," '{ if(($2 == 35) || ($2 == 80) || ($2)) { print } }' P1_output_file.csvi

rm -rf P1_filtered_*

for moduleID in 35 80 68 45 22 23 83 5 3 66 30 77 84; do
awk -F "," -v mod="$moduleID" '{ if(($2 == mod)) { print $1} }' P1_output_file.csv >> P1_filtered_TAM112A.csv
done

for moduleID in 28 20 38 46 52 25 90 4 65 53 91 55 37 88 72 11 1 29; do
awk -F "," -v mod="$moduleID" '{ if(($2 == mod)) { print $1} }' P1_output_file.csv >> P1_filtered_TAM112B.csv
done

for moduleID in 28 20 38 46 52 25 90 4 65 53 91 55 37 88; do
awk -F "," -v mod="$moduleID" '{ if(($2 == mod)) { print $1} }' P1_output_file.csv >> P1_filtered_TAM112B_I.csv
done

for moduleID in 4 65 53 91 55 37 88 72 11 1 29; do
awk -F "," -v mod="$moduleID" '{ if(($2 == mod)) { print $1} }' P1_output_file.csv >> P1_filtered_TAM112B_II.csv
done

for moduleID in 46 52 25 90 4 65 53 91 55 37 88 72 11 1 29; do
awk -F "," -v mod="$moduleID" '{ if(($2 == mod)) { print $1} }' P1_output_file.csv >> P1_filtered_TAM112B_III.csv
done