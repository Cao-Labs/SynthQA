# -*- mode: python -*-

block_cipher = None


a = Analysis(['assess_protein.py'],
             pathex=['protein_scoring'],
             binaries=None,
             datas=[],
             hiddenimports=['scipy._lib.messagestream', 'sklearn.neighbors.typedefs'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['matplotlib', 'PIL'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='sbrod',
          debug=False,
          strip=False,
          upx=True,
          console=True )
