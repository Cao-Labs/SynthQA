#!/bin/bash

# Usage: ./test_requests.sh localhost:8080

URL="$1"

curl -X POST \
-F target=T9999 \
-F sequence="VSYSDGHFLTKSGGVINFRKTRVTSITITILGNYGLRVVNGELQNTPLTFKGADFKSSTLKDELLIPLEGAVQLNTAPSTALCIFITTDHVYRELCMMQFLTDVDKTPFLVVLRSESKHETIQYMHIVTVHPFLSLT" \
-F tarball="http://predictioncenter.org/download_area/CASP12/server_predictions/T0860.stage1.3D.srv.tar.gz" \
-F email="xxxxxxx@gmail.com" \
-F model="plus" \
"$URL"
