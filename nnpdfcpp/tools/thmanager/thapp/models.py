 # This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin sqlcustom [app_label]'
# into your database.
from __future__ import unicode_literals

from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator

class Theory(models.Model):
    PTOS = (
        (0, 'LO'),
        (1, 'NLO'),
        (2, 'NNLO'),
        )
    FNSS = (
        ('FFNS','FFNS'),
        ('ZM-VFNS','ZM-VFNS'),
        ('FONLL-A','FONLL-A'),
        ('FONLL-B','FONLL-B'),
        ('FONLL-C','FONLL-C'),
        )
    MODEVS = (
        ('EXA','EXA'),
        ('EXP','EXP'),
        ('TRN','TRN'),
        )
    NFFFS = (
        (3,'3'),
        (4,'4'),
        (5,'5'),
        (6,'6'),
        )
    SXORDS = (
        ('LL','LL'),
        ('NLL','NLL'),
        ('NNLL','NNLL'),
        )
    HQS = (
        ('POLE','POLE'),
        ('MSBAR','MSBAR'),
        )
    
    id = models.IntegerField(db_column='ID', primary_key=True)
    pto = models.IntegerField(db_column='PTO', choices=PTOS, default=0) 
    fns = models.CharField(db_column='FNS', choices=FNSS, max_length=10, default='0')
    damp = models.BooleanField(db_column='DAMP', default=False)
    ic = models.BooleanField(db_column='IC', default=False)
    modev = models.CharField(db_column='ModEv', choices=MODEVS, max_length=3, default='TRN')
    xir = models.DecimalField(db_column='XIR', default=1.0, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)])  
    xif = models.DecimalField(db_column='XIF', default=1.0, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)])
    nfff = models.IntegerField(db_column='NfFF', choices=NFFFS, default='5')  
    maxnfas = models.IntegerField(db_column='MaxNfAs', choices=NFFFS, default='5') 
    maxnfpdf = models.IntegerField(db_column='MaxNfPdf', choices=NFFFS, default='5')  
    q0 = models.DecimalField(db_column='Q0', default=1.0, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)] )
    alphas = models.DecimalField(default=0.118, max_digits=5, decimal_places=4, validators=[MinValueValidator(0)])
    qref = models.DecimalField(db_column='Qref', default=91.2, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)] )
    qed = models.BooleanField(db_column='QED', default=False) 
    alphaqed = models.DecimalField(default=0.007496252,max_digits=15, decimal_places=9, validators=[MinValueValidator(0)])
    qedref = models.DecimalField(db_column='Qedref', default=1.777, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)] ) 
    sxres = models.BooleanField(db_column='SxRes', default=False) 
    sxord = models.CharField(db_column='SxOrd', choices=SXORDS, max_length=3, default='0')
    hq = models.CharField(db_column='HQ', choices=HQS, max_length=5, default='0')

    mc = models.DecimalField(default=1.275, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)])
    qmc = models.DecimalField(db_column='Qmc', default=1.275, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)]) 
    mb = models.DecimalField(default=4.18, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)])
    qmb = models.DecimalField(db_column='Qmb',default=4.18, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)])  
    mt = models.DecimalField(default=173.07, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)])
    qmt = models.DecimalField(db_column='Qmt',default=173.07, max_digits=10, decimal_places=3, validators=[MinValueValidator(0)]) 
    ckm = models.TextField(db_column='CKM', default='0.97428 0.22530 0.003470 0.22520 0.97345 0.041000 0.00862 0.04030 0.999152')  
    mz = models.DecimalField(db_column='MZ',default=91.1876, max_digits=10, decimal_places=4, validators=[MinValueValidator(0)]) 
    mw = models.DecimalField(db_column='MW',default=80.398, max_digits=10, decimal_places=4, validators=[MinValueValidator(0)])
    gf = models.DecimalField(db_column='GF', default=1.1663787e-05, max_digits=15, decimal_places=15, validators=[MinValueValidator(0)])  # Field name made lowercase.
    sin2tw = models.DecimalField(db_column='SIN2TW', default=0.23126, max_digits=15, decimal_places=10, validators=[MinValueValidator(0)])  # Field name made lowercase.
    tmc = models.BooleanField(db_column='TMC', default=True)  # Field name made lowercase.
    mp = models.DecimalField(db_column='MP', default=0.938, max_digits=5, decimal_places=3, validators=[MinValueValidator(0)])  # Field name made lowercase.

    comments = models.TextField(db_column='Comments', blank=False, null=False)

    class Meta:
        managed = False
        db_table = 'TheoryIndex'
