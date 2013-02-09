#!/bin/bash

echo "Polyhedra performance test with car example."
echo "POLYHEDRA PERFORMANCE TEST WITH CAR EXAMPLE." > log
for i in {2..50}
do
echo "Do for $i states."
echo "DO FOR $i STATES" >> log
python car.py $i
./../../prism/bin/prism car.gen.smg car.gen.props >> log
done
echo "Done, exiting ..."
