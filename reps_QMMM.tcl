proc rep_charged {} {
    mol delrep 0 top
    mol representation CPK 1.00000 0.300000 40.000000 40.000000
    mol color element
    mol selection "index 1634 1635 1636 1637 1638 1639 1640 1641 1642 1643 
    1644 1645 1646 1647 1648 1649 1650 1651 1652 1653 1654 1655 1656 1657 
    1658 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671 
    5852 5853 5854 5855 5856 5857 5858 5859 5860 5861 5862 5863 5864 7126 
    7127 7128 7129 7130 7131 7132 7133 7134 7135 7136 7137 7138 9960 9961 
    9962 9963 9964 9965 9966 9967 9968 9969 9970 9971 9972 10337 10338 10339 
    10340 10341 10342 10343 10344 10345 10346 10347 10348 10349 12651 12652 12653 
    12654 12655 12656 12657 12658 12659 12660 12661 12662 12663 30773 30774 30775 
    30776 30777 30778 30779 30780 30781 30782 30783 30784 30785 38625 38626 38627 
    38628 38629 38630 38631 38632 38633 38634 38635 38636 38637 59670 59671 59672 
    59673 59674 59675 59676 59677 59678 59679 59680 59681 59682 59683 64318 64319 
    64320 64321 64322 64323 64324 64325 64326 64327 64328 64329 64330 64331 64332 
    64333 64334 64335 64336 64337 64338 68541 70214 70215 70216 70217 70218 70219 
    70220 70221 70222 70223 70224 70225 70226 70227 70228 70229 70230 70231 70232 
    70233 70234 78158 78159 78160 78161 78162 78163 78164 78165 78166 78167 78168 
    78169 78170 78171 83240 83241 83242 83243 83244 83245 83246 83247 83248 83249 
    83250 83251 83252 83253 83254 83255 83256"
    mol addrep top
    mol representation CPK 1.00000 0.300000 40.000000 40.000000
    mol selection "(same resid as within 6 of index 1636) and not index 1634 1635 1636 1637 1638 1639 1640 1641 1642 1643 
    1644 1645 1646 1647 1648 1649 1650 1651 1652 1653 1654 1655 1656 1657 
    1658 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671 
    5852 5853 5854 5855 5856 5857 5858 5859 5860 5861 5862 5863 5864 7126 
    7127 7128 7129 7130 7131 7132 7133 7134 7135 7136 7137 7138 9960 9961 
    9962 9963 9964 9965 9966 9967 9968 9969 9970 9971 9972 10337 10338 10339 
    10340 10341 10342 10343 10344 10345 10346 10347 10348 10349 12651 12652 12653 
    12654 12655 12656 12657 12658 12659 12660 12661 12662 12663 30773 30774 30775 
    30776 30777 30778 30779 30780 30781 30782 30783 30784 30785 38625 38626 38627 
    38628 38629 38630 38631 38632 38633 38634 38635 38636 38637 59670 59671 59672 
    59673 59674 59675 59676 59677 59678 59679 59680 59681 59682 59683 64318 64319 
    64320 64321 64322 64323 64324 64325 64326 64327 64328 64329 64330 64331 64332 
    64333 64334 64335 64336 64337 64338 68541 70214 70215 70216 70217 70218 70219 
    70220 70221 70222 70223 70224 70225 70226 70227 70228 70229 70230 70231 70232 
    70233 70234 78158 78159 78160 78161 78162 78163 78164 78165 78166 78167 78168 
    78169 78170 78171 83240 83241 83242 83243 83244 83245 83246 83247 83248 83249 
    83250 83251 83252 83253 83254 83255 83256"
    mol color colorid 0
    mol addrep top
    mol representation CPK 1.00000 0.300000 40.000000 40.000000
    mol selection "(same resid as index 1634 1635 1636 1637 1638 1639 1640 1641 1642 1643 
    1644 1645 1646 1647 1648 1649 1650 1651 1652 1653 1654 1655 1656 1657 
    1658 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671 
    5852 5853 5854 5855 5856 5857 5858 5859 5860 5861 5862 5863 5864 7126 
    7127 7128 7129 7130 7131 7132 7133 7134 7135 7136 7137 7138 9960 9961 
    9962 9963 9964 9965 9966 9967 9968 9969 9970 9971 9972 10337 10338 10339 
    10340 10341 10342 10343 10344 10345 10346 10347 10348 10349 12651 12652 12653 
    12654 12655 12656 12657 12658 12659 12660 12661 12662 12663 30773 30774 30775 
    30776 30777 30778 30779 30780 30781 30782 30783 30784 30785 38625 38626 38627 
    38628 38629 38630 38631 38632 38633 38634 38635 38636 38637 59670 59671 59672 
    59673 59674 59675 59676 59677 59678 59679 59680 59681 59682 59683 64318 64319 
    64320 64321 64322 64323 64324 64325 64326 64327 64328 64329 64330 64331 64332 
    64333 64334 64335 64336 64337 64338 68541 70214 70215 70216 70217 70218 70219 
    70220 70221 70222 70223 70224 70225 70226 70227 70228 70229 70230 70231 70232 
    70233 70234 78158 78159 78160 78161 78162 78163 78164 78165 78166 78167 78168 
    78169 78170 78171 83240 83241 83242 83243 83244 83245 83246 83247 83248 83249 
    83250 83251 83252 83253 83254 83255 83256) and not within 6 of index 1636"
    mol color colorid 1
    mol addrep top

}

proc rep_charged_8 {} {
    mol delrep 0 top
    mol representation CPK 1.00000 0.300000 40.000000 40.000000
    mol color element
    mol selection "index 1634 to 1671 5852 to 5864 7126 to 7138 9960 to 9972 
    10337 to 10349 12313 to 12325 12404 to 12416 12651 to 12663 13210 to 13222 
    19931 to 19943 28446 to 28458 30773 to 30785 32710 to 32722 33958 to 33970 
    38625 to 38637 41264 to 41276 43149 to 43161 59670 to 59683 62344 to 62357 
    64318 to 64338 64779 to 64779 68541 to 68541 70214 to 70234 78158 to 78171 
    79712 to 79725 83240 to 83256 "
    mol addrep top
}