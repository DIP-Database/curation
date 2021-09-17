
filename='pdbListHeader.txt'
exec 4<"$filename"

while read -u4 p ; do
    echo Beginning: "$p" >> log.txt
    echo Beginning: "$p"
    python3 pdbAnalysis.py -s "$p" >> log.txt
    echo Completed: "$p" >> log.txt
    echo Completed: "$p"
done
