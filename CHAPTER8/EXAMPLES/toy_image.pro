PRO toy_image
 image = fltarr (800 ,800 ,3)
 b = 2.0
 image [99:699 ,299:499 ,*] = b
 image [299:499 ,99:699 ,*] = b
 image [299:499 ,299:499 ,*] = 2*b
 image = smooth (image ,[13 ,13 ,1] ,/ edge_truncate )
 n1 = randomu (seed ,800 ,800 ,/ normal )
 n2 = randomu (seed ,800 ,800 ,/ normal )
 n3 = randomu (seed ,800 ,800 ,/ normal )
 image [* ,* ,0] = image [* ,* ,0] + n1
 image [* ,* ,1] = image [* ,* ,1] + n2 + n1
 image [* ,* ,2] = image [* ,* ,2] + n3 + n1 /2 + n2 /2
 envi_enter_data , image
end