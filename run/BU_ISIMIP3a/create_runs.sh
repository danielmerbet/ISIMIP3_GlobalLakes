c=0
for p in {1..41449..415}; do
  let c++
  cp run.sh runs/run${p}.sh
  sed -i "26s/.*/for pixel in {${p}..$((p+414))}/" runs/run${p}.sh
  sed -i "21s/.*/folder=${c}/" runs/run${p}.sh
  
done
