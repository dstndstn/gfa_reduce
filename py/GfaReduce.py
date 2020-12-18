#!/usr/bin/env python

from DOSlib.application import Application
import DOSlib.discovery as discovery

class GfaReduce(Application):
    commands = ['configure', 'get', 'set', 'gfaproc']
    defaults = {'simulate': False,
                'gfa_reduce_etc_dir': "${GFA_REDUCE_META}/",
                }

    def init(self):
        self.loglevel('INFO')
        self.info('INITIALIZING Role %s' % self.role)

        self.instrument_dir = self.config['gfa_reduce_etc_dir']

        # Create status shared variable
        self.status_sv = self.shared_variable('STATUS')
        self.status_sv.publish(republish=True)
        self.status_sv.write('INITIALIZING')

        # OCS Messages
        self.OCSMessagesSV = self.shared_variable("MESSAGES", group="OCS")
        self.OCSMessagesSV.publish(allowMultiplePublishers=True)

        # use the discovery module to announce to the world
        if self.connected:
            self._setup_discovery(discovery.discoverable)
        self.status_sv.write('INITIALIZED')
        self.info('INITIALIZED')
        return self.SUCCESS

    def _setup_discovery(self, discoverable):
        self.info('setting up discovery')
        discoverable(role=self.role, tag='GFAREDUCE', interface=self.role)

    # PML functions (functions that can be accessed remotely)
    def configure(self, constants='DEFAULT'):
        # reset status
        if self.status_sv._value == 'INITIALIZING':
            raise RuntimeError("init function not completed yet, unable to configure")
        # reset application status
        if self.status_sv._value != 'READY':
            self.status_sv.write('CONFIGURING')
        # ...
        self.info('CONFIGURED')
        self.status_sv.write('READY')
        return self.SUCCESS

    def get(self, param, config_id=None):
        """Retrieve parameters dynamically.

        Options include:
        """
        # fix inconsistent variable names
        if param.lower() in self.getter_setter_handles:
            param = self.getter_setter_handles[param.lower()]

        if param == 'status':
            value = self.status_sv._value
        elif hasattr(self, param):
            value = getattr(self, param)
        elif param in self.config:
            value = self.config[param]
        else:
            return 'Invalid parameter for get command: %s' % param
        return value

    def set(self, *args, **params):
        """Set parameters dynamically."""
        # for param in params:
        #     value = params[param]
        # 
        #     # fix inconsistent variable names
        #     if param.lower() in self.getter_setter_handles:
        #         param = self.getter_setter_handles[param.lower()]
        # 
        #     if self.dervish_plate_maker \
        #        and hasattr(self.dervish_plate_maker, param) or param in handle_synonyms:
        #             setattr(self.dervish_plate_maker, param, value)
        #             if param in ['instrument', 'inst_name']:
        #                 self.inst_name = value
        #     elif param in ['instrument', 'inst_name']:
        #         self.inst_name = value
        #     elif hasattr(self, param):
        #         setattr(self, param, value)
        #     else:
        #         return 'FAILED: Invalid parameter for set command: %s' % params
        return self.SUCCESS

    def gfaproc(self, *args, **params):
        """Run gfa_reduce.

        Returns:
            If successful, the string "SUCCESS"
            otherwise, a string starting with "FAILED"
        """
        self.info('calling gfa_reduce')
        self.set(*args, **params)

        self.gfa_center = None
        self.gfa_guide  = None
        self.guider_wcs = None
        
        import gfa_reduce.gfa_red as gfa_red
        from collections import OrderedDict
        import numpy as np
        fm = gfa_red.acquire_field(gfa_targets=gfa_targets, exp_data=hdulist)

        #### FIXME -- hexrate??

        self.gfa_center = OrderedDict([('ra', fm.ra),
                                       ('dec', fm.dec),
                                       ('hexrot', fm.hexrot_deg),
                                       ('hexrate', 0.),])
        #self.gfa_guide = ...guide stars...

        wcsvals = []
        for hdu in hdulist:
            hdr = hdu.header
            # Strip "GUIDE" off "GUIDE2" to return 2.
            ext = hdr['EXTNAME']
            if not ext.startswith('GUIDE'):
                continue
            guide_loc = int(ext[5:], 10)
            wcsvals.append(tuple([guide_loc] + [hdr[k] for k in [
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']]))
        self.guider_wcs = np.rec.array(wcsvals,
             dtype=np.dtype([('GFA_LOC', '<i8'), ('CRVAL1', '<f8'), ('CRVAL2', '<f8'),
                             ('CRPIX1', '<f8'), ('CRPIX2', '<f8'),
                             ('CD1_1', '<f8'), ('CD1_2', '<f8'),
                             ('CD2_1', '<f8'), ('CD2_2', '<f8')]))

        if result.startswith('FAILED'):
            self.error(result)
            return result
        return self.SUCCESS

    # Instance Connection Callbacks
    def about_to_connect_to_instance(self, *args, **kwargs):
        pass

    def did_connect_to_instance(self, *args, **kwargs):
        self.info('connected, setting up discovery stuff')
        discovery.reset()
        discovery.reset_discovered()
        self._setup_discovery(discovery.discoverable)

    def about_to_disconnect_from_instance(self, *args, **kwargs):
        pass

    def did_disconnect_from_instance(self, *args, **kwargs):
        self.info('disconnected, clearing discovery stuff')
        discovery.reset()
        discovery.reset_discovered()

    def main(self):
        while not self.shutdown_event.is_set():
            time.sleep(1)
        print('GfaReduce exits')

if __name__ == '__main__':
    GfaReduce().run()
