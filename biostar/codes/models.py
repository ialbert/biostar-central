from django.db import models

MAX_TEXT_FIELD = 256
MAX_FLOAT_FIELD = 256
MAX_DECIMALS_DIGITS = 20


class Location(models.Model):

    name = models.CharField(max_length=MAX_TEXT_FIELD)

    # Unique site identifier
    site_id = models.CharField(max_length=MAX_TEXT_FIELD, unique=True)

    latitude = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)
    logitude = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)

    def __str__(self):
        return self.name


class QPCRSample(models.Model):
    """Data from sample sheet"""

    # Site where sample is collected
    location = models.ForeignKey(Location, null=True, on_delete=models.SET_NULL)

    # Unique sample identifier
    sample_id = models.CharField(max_length=MAX_TEXT_FIELD, unique=True)

    type = models.CharField(max_length=MAX_TEXT_FIELD)

    #TODO: diffrence between this and extraction_date?
    date = models.DateField(null=True)
    time = models.TimeField(null=True)

    # Total number of filters
    nfilters = models.IntegerField(null=True)

    # Comma separated filter volumes
    filters = models.CharField(null=True, max_length=MAX_TEXT_FIELD)

    total_volume = models.IntegerField(null=True)
    extraction_date = models.DateField(null=True)

    elution_volume = models.IntegerField(null=True)
    pcr_template_volume = models.IntegerField(null=True)
    pcr_concentration = models.IntegerField(null=True)

    # Water temperature
    water_temp = models.FloatField(max_length=MAX_FLOAT_FIELD, null=True)

    def __str__(self):
        return self.sample_id


class QPCRExperiment(models.Model):

    name = models.CharField(max_length=MAX_TEXT_FIELD)

    barcode = models.CharField(max_length=MAX_TEXT_FIELD, unique=True)

    target_name = models.CharField(max_length=MAX_TEXT_FIELD, null=True)

    # All samples present in the experiment
    samples = models.ManyToManyField(QPCRSample, related_name="samples")

    instrument_serial_number = models.IntegerField()
    instrument_type = models.CharField(max_length=MAX_TEXT_FIELD)

    user = models.CharField(max_length=MAX_TEXT_FIELD)

    def add_samples(self, samples):
        """Add more samples to experiment """

        # Remove samples then add to ensure no duplicates
        self.samples.remove(*samples)
        self.samples.add(*samples)

    def __str__(self):
        return self.name


class QPCRMeasurement(models.Model):
    """Models qPCR outputs """

    # The well position (A1, B3, etc. )
    position = models.CharField(max_length=MAX_TEXT_FIELD)

    # Particular sample found in this well position.
    sample = models.ForeignKey(QPCRSample, on_delete=models.CASCADE)

    # Experiment plate this well belongs to.
    experiment = models.ForeignKey(QPCRExperiment, on_delete=models.CASCADE)

    ct = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)
    ct_mean = models.FloatField(max_length=MAX_FLOAT_FIELD, null=True)
    ct_sd = models.FloatField(max_length=MAX_FLOAT_FIELD, null=True)

    quantity = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)
    quantity_mean = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)
    quantity_sd = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)

    automatic_ct_threshold = models.BooleanField(default=False)
    ct_threshold = models.DecimalField(max_digits=MAX_DECIMALS_DIGITS, null=True, decimal_places=12)

    automatic_baseline = models.BooleanField(default=False)
    baseline_start = models.IntegerField(null=True)
    baseline_end = models.IntegerField(null=True)

    def __str__(self):
        return f"{self.position}-{self.experiment}"