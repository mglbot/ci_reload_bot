#!/bin/bash

echo ""
echo ""
echo "checking the code against clang-format"
echo ""
echo ""

cd src/

STATUS=0

for bob in $(find . -name "*.cpp" -o -name "*.hpp" -o -name "*.c" -o -name "*.h")
do
  DIFF=$(diff -u <(cat ${bob}) <(/usr/bin/clang-format-6.0 -style=file ${bob}))
  if [ ! -z "${DIFF}" ] ; then
    STATUS=1
    echo ${bob}
    echo "${DIFF}"
  fi
done

exit ${STATUS}
