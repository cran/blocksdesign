# ************************************************************************************************************************
# Orthogonal basis for M internal function only
# ************************************************************************************************************************
orthogM=function(M) {
  QR = qr(M) # qr transformation
  M = qr.Q(QR)[,1:QR$rank,drop=FALSE] # orthogonal basis for M where Q'Q=I
  M
}
