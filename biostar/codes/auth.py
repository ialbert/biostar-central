import os
import csv
import datetime
import zipfile
from functools import partial

# Try import models without setting up django
from biostar.codes.models import QPCRSample, QPCRMeasurement, Location, QPCRExperiment


class Bunch(object):
    def __init__(self, **kwargs):

        self.site_id = self.site_description = self.longitude = self.latitude = None
        self.date = self.time = self.extraction_date = self.sample_id = None

        self.sample_type = self.total_volume = self.nfilters = self.elution_volume = None
        self.pcr_template_volume = self.pcr_concentration = self.water_temp = None

        self.position = self.ct = self.ct_mean =self.ct_sd = self.quantity = None
        self.quantity_mean = self.quantity_sd = self.ct_threshold =None

        self.baseline_start = self.baseline_end = self.target_name = self.sample_id = None
        self.automatic_ct_threshold = self.automatic_baseline =None

        self.__dict__.update(**kwargs)


def create_location(site_id, site_name, longitude, latitude):
    """Creates a Location object in the database if it does not exist."""

    location = Location.objects.filter(site_id=site_id).first()
    if location is None:
        location = Location.objects.create(site_id=site_id, name=site_name,
                                           latitude=latitude, logitude=longitude)
    return location


def create_qpcr_sample(location, sample_id, sample_type, total_volume, date, time,
                       nfilters, extraction_date, elution_volume, pcr_template_volume,
                       pcr_concentration, water_temp):
    """Creates a Sample object in the database if it does not exist."""

    sample = QPCRSample.objects.filter(sample_id=sample_id).first()
    if sample is None:
        sample = QPCRSample.objects.create(location=location, sample_id=sample_id,
                                           type=sample_type, total_volume=total_volume,
                                           date=date, time=time, nfilters=nfilters,
                                           extraction_date=extraction_date, elution_volume=elution_volume,
                                           pcr_template_volume=pcr_template_volume,
                                           pcr_concentration=pcr_concentration,
                                           water_temp=water_temp)
    return sample


def create_qpcr_experiment(**info):
    """Creates an experiment.
    Returns a tuple with the object and True or False stating if it already existed.

    """

    barcode = info["experiment_barcode"]
    exp = QPCRExperiment.objects.filter(barcode=barcode).first()
    already_exists = True

    if exp is None:
        name = info["experiment_name"]
        user = info["experiment_user_name"]
        serial_number = info["instrument_serial_number"]
        instrument_type = info["instrument_type"]
        exp = QPCRExperiment.objects.create(barcode=barcode, name=name, user=user,
                                            instrument_serial_number=serial_number,
                                            instrument_type=instrument_type)
        already_exists = False

    return exp, already_exists


def create_qpcr_measurement(experiment, sample, position, ct, ct_mean, ct_sd, quantity_mean,
                            quantity, quantity_sd, automatic_ct_threshold, ct_threshold, automatic_baseline,
                            baseline_start, baseline_end):
    """Creates measurement """

    meas = QPCRMeasurement.objects.create(experiment=experiment, sample=sample,
                                          position=position, ct=ct, ct_mean=ct_mean,
                                          ct_sd=ct_sd, quantity_mean=quantity_mean,
                                          quantity=quantity, quantity_sd=quantity_sd,
                                          automatic_ct_threshold=automatic_ct_threshold, ct_threshold=ct_threshold,
                                          automatic_baseline=automatic_baseline, baseline_start=baseline_start,
                                          baseline_end=baseline_end)
    return meas


def fill_header(first_row, header={}):
    """Fill header dict with {column name: column index} """

    for idx in range(len(first_row)):
        # clean the column key
        key = filter(lambda x: x.isalnum(), first_row[idx].split())
        key = "_".join(list(key)).lower()
        header[key] = idx

    return header


def parse_date(date):

    if date:
        pattern = "%m/%d/%Y"
        return datetime.datetime.strptime(date, pattern).date()


def bunch_samplesheet_row(row, header):
    """Create a bunch object for a single row using header of sample sheet"""

    parse_time = lambda time: None if time is None else time.split()[0].strip()

    bunched = Bunch(
        site_id=row[header["site_id"]],
        site_description=row[header["site_description"]],
        longitude=row[header["longitude"]],
        latitude=row[header["latitude"]],
        date=parse_date(row[header["date"]]),
        time=parse_time(row[header["time"]]),
        extraction_date=parse_date(row[header["extraction_date"]]),
        sample_id=row[header["sample_id"]],
        sample_type=row[header["sample_type"]],
        total_volume=row[header["total_volume"]],
        nfilters=row[header["total_number_of_filters"]],
        elution_volume=row[header["elution_volume"]],
        pcr_template_volume=row[header["pcr_template_volume"]],
        pcr_concentration=row[header["pcr_concentration"]],
        water_temp=row[header["water_temp"]]
          )

    return bunched


def load_samplesheet(sample_sheet):
    """Create Location and Sample objects from csv file."""

    csv_reader = csv.reader(open(sample_sheet, "r"), delimiter=",")
    header = fill_header(first_row=next(csv_reader))

    for row in csv_reader:

        # Change NA values to None
        row = list(map(lambda val: None if val == "NA" else val, row))
        row = bunch_samplesheet_row(row=row, header=header)

        # Enter location
        location = create_location(site_id=row.site_id,site_name=row.site_description,
                                   longitude=row.longitude, latitude=row.latitude)
        # Enter meta-data
        create_qpcr_sample(location=location, sample_id=row.sample_id, sample_type=row.sample_type,
                           total_volume=row.total_volume, date=row.date, time=row.time,
                           nfilters=row.nfilters, extraction_date=row.extraction_date,
                           elution_volume=row.elution_volume, pcr_template_volume=row.pcr_template_volume,
                           pcr_concentration=row.pcr_concentration, water_temp=row.water_temp
                           )


def parse_exp_info(stream):
    """Parse an experiment info and header."""

    exp_info = dict()
    header = dict()
    for row in stream:
        if row[0].startswith("Well"):
            header = fill_header(first_row=row, header=header)
            break
        # Clean the key to avoid mistakes when looking up later
        key = filter(lambda x: x.isalnum(), row[0].split())
        key = "_".join(list(key)).lower()
        exp_info[key] = row[1]

    return exp_info, header


def bunch_measurement_row(row, header):

    """Create a bunch object for a single row using header from qPCR output"""

    bool_map = dict(true=True, false=False, Y=True, N=False)

    bunched = Bunch(
        position=row[header["well_position"]],
        ct=row[header["ct"]],
        ct_mean=row[header["ct_mean"]],
        ct_sd=row[header["ct_sd"]],
        quantity=row[header["quantity"]],
        quantity_mean=row[header["quantity_mean"]],
        quantity_sd=row[header["quantity_sd"]],
        ct_threshold=row[header["ct_threshold"]],
        baseline_start=row[header["baseline_start"]],
        baseline_end=row[header["baseline_end"]],
        target_name=row[header["target_name"]],
        sample_id=row[header["sample_name"]],

        # Boolean values
        automatic_ct_threshold=bool_map[row[header["automatic_ct_threshold"]]],
        automatic_baseline=bool_map[row[header["automatic_baseline"]]],

    )

    return bunched


def load_experiment(csvfile):
    """Creates Experiment and Measurement objects from the qPCR measurement."""

    measurements = csv.reader(open(csvfile, "r"), delimiter=",")
    exp_info, header = parse_exp_info(stream=measurements)
    experiment, already_exists = create_qpcr_experiment(**exp_info)
    if already_exists:
        # Does not read further into file if experiment already exists.
        return

    for row in measurements:

        # Clean row of empty values
        is_empty = lambda x: x == "UNKNOWN" or len(x) == 0 or x == "Undetermined"
        row = list(map(lambda val: None if is_empty(val) else val.replace(",", ""), row))
        row = bunch_measurement_row(row=row, header=header)
        sample = QPCRSample.objects.filter(sample_id=row.sample_id).first()

        # Add sample to experiment
        experiment.add_samples(samples=[sample])

        # Set the target name if it has not been set yet
        if not experiment.target_name:
            experiment.target_name = row.target_name

        # Enter the measurement data
        create_qpcr_measurement(experiment=experiment, sample=sample,
                                position=row.position, ct=row.ct, ct_mean=row.ct_mean,
                                ct_sd=row.ct_sd, quantity_mean=row.quantity_mean,
                                quantity=row.quantity, quantity_sd=row.quantity_sd,
                                automatic_ct_threshold=row.automatic_ct_threshold,
                                ct_threshold=row.ct_threshold, automatic_baseline=row.automatic_baseline,
                                baseline_start=row.baseline_start, baseline_end=row.baseline_end)

    # Commit changes made
    experiment.save()


def load_all(sample_sheet, data_dir=None, zip_file=None):
    """Load sample sheet and data into the database."""

    # Load the sample independent of other things.
    load_samplesheet(sample_sheet=sample_sheet)

    # Iterate over files and load into database
    results = [] if not data_dir else os.scandir(data_dir)

    # Load experiments
    list(map(load_experiment, results))

    # Unpack and load zipfile
    if zip_file:

        zip_ref = zipfile.ZipFile(zip_file, 'r')
        files = list(filter(lambda x: not x.filename.startswith("__"), zip_ref.infolist()))

        extract = partial(zip_ref.extract, path=os.path.join(BASE_DIR, "data"))

        # Extract and return a list of files.
        output = list(filter(lambda f: os.path.isfile(f), map(extract, files)))

        list(map(load_experiment, output))
